function [F, G, dg] = cpmg_dynamics(x, ug, d, p, tie_in)
%CPMG_DYNAMICS  Nonlinear control-affine dynamics of one cyber-physical
%   microgrid area.
%
%   [F, G, dg] = CPMG_DYNAMICS(x, ug, d, p, tie_in) returns the drift field
%   F and the edge input matrix G such that
%
%       xdot = F(x, ug, d) + G(x) * u_edge ,        u_edge = [u_B ; u_A]
%
%   The system is control-affine IN THE EDGE INPUT ONLY.  The central AGC
%   command ug is treated as an exogenous, zero-order-held signal and is
%   folded into F, because during a DoS event it is frozen and is therefore
%   part of the known drift as seen by the edge controller.
%
%   INPUTS
%     x      9x1  canonical state, indexed by cpmg_spec().ix
%     ug     1x1  central AGC command to the governor reference   [p.u.]
%     d      3x1  [dP_res ; dP_L ; dT_out]
%     p      1x1  parameter struct from params_isolated / params_39bus
%     tie_in 1x1  tie-line coupling term 2*pi*sum_j T_ij*(df_i - df_j).
%                 Pass 0 (or omit) for an isolated microgrid.
%
%   OUTPUTS
%     F      9x1  drift vector field
%     G      9x2  edge input matrix, CONTAINS Phi(SoC) in the BESS column
%     dg     struct of diagnostics required by the barrier module and by
%            verify_plant.m:
%              .Phi        SoC-dependent BESS derating in [0,1]
%              .Phi_degen  true when Phi < p.Phi_min (relative degree loss)
%              .PL         composite load power                    [p.u.]
%              .dPL_ddf    d(PL)/d(df), needed for Lf^2 h          [p.u./p.u.]
%              .dPm_dot    actual (GRC- and limit-clipped) dPm/dt  [p.u./s]
%              .ACE        beta*df + dPtie                         [p.u.]
%              .db_active  governor dead-band is suppressing droop
%              .grc_active GRC rate limit is clipping
%              .lim_active [Pv at limit , Pm at limit]
%              .COP_eff    effective coefficient of performance
%
%   MODEL, in the sign convention frozen by cpmg_spec (D3)
%   ------------------------------------------------------
%   (1) Swing, nonlinear composite load:
%         M df' = dPm + dPB + dPres - dPA - dPtie - [P_L(df,dPL) - PL0]
%         P_L(df,dPL) = (PL0 + dPL)*(1 + D*df)^gamma
%
%   (2) Tie line (linear synchronising-power model, Assumption 3):
%         dPtie' = tie_in
%
%   (3) Turbine with generation rate constraint:
%         dPm' = sat_GRC( (-dPm + dPv)/Tt )
%
%   (4) Governor with dead-band:
%         dPv' = ( -dPv - phi_db(df)/R + ug )/Tg
%
%   (5) Area control error integral:
%         zACE' = beta*df + dPtie
%
%   (6) BESS with SoC-dependent saturation (D4):
%         dPB' = -dPB/TB + Phi(SoC)*u_B/TB
%
%   (7) State of charge:
%         SoC' = -eta_B*dPB/E_B          (discharge positive -> SoC falls)
%
%   (8) AIAC aggregate consumption:
%         dPA' = -dPA/Tac + u_A/Tac
%
%   (9) Aggregate building temperature (ETP, after Hui2019):
%         dTin' = a_T*(Tout - Tin) - b_T*(PA_nom + dPA)
%
%   RELATIVE DEGREE NOTE
%   --------------------
%   With h_f(x) = df - df_min, the AIAC and BESS inputs first appear at the
%   SECOND derivative of h_f, giving uniform relative degree 2 provided
%   Phi > 0.  The governor dead-band acts on dPv, which appears only at the
%   THIRD derivative, so the HO-CBF construction never differentiates through
%   the dead-band and no smooth surrogate is required for control.  The GRC
%   does appear at the second derivative (through dPm_dot) and is therefore
%   returned exactly, clipping included, in dg.dPm_dot.
%
%   See also CPMG_SPEC, PARAMS_ISOLATED, PARAMS_39BUS, CPMG_DYNAMICS_MULTIAREA.

if nargin < 5 || isempty(tie_in), tie_in = 0; end

spec = p.spec;
ix   = spec.ix;
eps0 = spec.num.eps_guard;

%% ------------------------------------------------------------------ %%
%  Unpack
%% ------------------------------------------------------------------ %%
df    = x(ix.df);
dPtie = x(ix.dPtie);
dPm   = x(ix.dPm);
dPv   = x(ix.dPv);
dPB   = x(ix.dPB);
SoC   = x(ix.dSB);
dPA   = x(ix.dPA);
dTin  = x(ix.dTin);

dPres = d(spec.id.dPres);
dPL   = d(spec.id.dPL);
dTout = d(spec.id.dTout);

F = zeros(spec.nx,1);
G = zeros(spec.nx,spec.nu);

%% ------------------------------------------------------------------ %%
%  (1) Composite frequency-dependent load  P_L(df, dPL)
%% ------------------------------------------------------------------ %%
if spec.nl.load_frequency
    base    = max(1 + p.D*df, eps0);                     % guard: no complex powers
    PL      = (p.PL0 + dPL)*base^p.gamma;
    dPL_ddf = (p.PL0 + dPL)*p.gamma*p.D*base^(p.gamma-1);
else
    % Linear fallback, retained only for the ablation study.  Never used in
    % the reported results.
    PL      = (p.PL0 + dPL) + p.D*df;
    dPL_ddf = p.D;
end

%% ------------------------------------------------------------------ %%
%  Swing equation
%% ------------------------------------------------------------------ %%
F(ix.df) = ( dPm + dPB + dPres - dPA - dPtie - (PL - p.PL0) )/p.M;

%% ------------------------------------------------------------------ %%
%  (2) Tie line
%% ------------------------------------------------------------------ %%
F(ix.dPtie) = tie_in;

%% ------------------------------------------------------------------ %%
%  (3) Turbine with generation rate constraint
%% ------------------------------------------------------------------ %%
raw_Pm = (-dPm + dPv)/p.Tt;
if spec.nl.turbine_GRC
    switch lower(spec.nl.grc_mode)
        case 'hard'
            dPm_dot = min(max(raw_Pm, -p.GRC), p.GRC);
        case 'smooth'
            dPm_dot = p.GRC*tanh(p.grc_smooth_k*raw_Pm/p.GRC);
        otherwise
            error('cpmg_dynamics:grcMode','Unknown grc_mode "%s"', spec.nl.grc_mode);
    end
else
    dPm_dot = raw_Pm;
end
grc_active = abs(raw_Pm) > p.GRC*(1 - 1e-9) && spec.nl.turbine_GRC;

% Anti-windup on the mechanical power magnitude limit
[dPm_dot, sig_Pm] = clampRate(dPm, dPm_dot, p.Pm_min, p.Pm_max, p.aw_band);
F(ix.dPm) = dPm_dot;

%% ------------------------------------------------------------------ %%
%  (4) Governor with dead-band
%% ------------------------------------------------------------------ %%
if spec.nl.governor_deadband
    switch lower(spec.nl.deadband_mode)
        case 'hard'
            phi_db = deadbandHard(df, p.df_db);
        case 'smooth'
            phi_db = deadbandSmooth(df, p.df_db, p.db_smooth_k);
        otherwise
            error('cpmg_dynamics:dbMode','Unknown deadband_mode "%s"', spec.nl.deadband_mode);
    end
else
    phi_db = df;
end
db_active = spec.nl.governor_deadband && abs(df) <= p.df_db;

raw_Pv = ( -dPv - phi_db/p.R + ug )/p.Tg;
[dPv_dot, sig_Pv] = clampRate(dPv, raw_Pv, p.Pv_min, p.Pv_max, p.aw_band);
F(ix.dPv) = dPv_dot;

%% ------------------------------------------------------------------ %%
%  (5) Area control error integral
%% ------------------------------------------------------------------ %%
ACE = p.beta*df + dPtie;
F(ix.zACE) = ACE;

%% ------------------------------------------------------------------ %%
%  (6)-(7) BESS with SoC-dependent saturation Phi   [D4]
%% ------------------------------------------------------------------ %%
if spec.nl.bess_Phi && spec.use_Phi_in_G
    Phi = socSaturation(SoC, p.SoC_lo, p.SoC_hi, p.ks);
else
    Phi = 1.0;
end
Phi_degen = Phi < p.Phi_min;

F(ix.dPB)      = -dPB/p.TB;
G(ix.dPB, spec.iu.uB) = Phi/p.TB;            % <-- Phi appears in LgLf_h

F(ix.dSB) = -p.eta_B*dPB/p.E_B_s;

% Hard energy bound: SoC cannot leave [0,1] regardless of Phi.  With the
% correct Phi this clamp should never fire; it is kept as an assertion trap.
if (SoC >= 1 && F(ix.dSB) > 0) || (SoC <= 0 && F(ix.dSB) < 0)
    F(ix.dSB) = 0;
end

%% ------------------------------------------------------------------ %%
%  (8) AIAC aggregate consumption
%% ------------------------------------------------------------------ %%
F(ix.dPA)      = -dPA/p.Tac;
G(ix.dPA, spec.iu.uA) = 1/p.Tac;

%% ------------------------------------------------------------------ %%
%  (9) Aggregate building thermal dynamics
%% ------------------------------------------------------------------ %%
Tin_abs  = p.Tin_ref  + dTin;
Tout_abs = p.Tout_nom + dTout;

if spec.nl.cop_temperature
    COP_eff = p.COP*(1 - p.COP_slope*((Tout_abs - Tin_abs) - p.COP_dTnom));
    COP_eff = min(max(COP_eff, p.COP_min), p.COP_max);
else
    COP_eff = p.COP;
end
bT_eff = p.b_T*(COP_eff/p.COP);

F(ix.dTin) = p.a_T*(Tout_abs - Tin_abs) - bT_eff*(p.PA_nom + dPA);

%% ------------------------------------------------------------------ %%
%  Diagnostics
%% ------------------------------------------------------------------ %%
dg = struct( ...
    'Phi',        Phi, ...
    'Phi_degen',  Phi_degen, ...
    'PL',         PL, ...
    'dPL_ddf',    dPL_ddf, ...
    'dPm_dot',    dPm_dot, ...
    'dPv_dot',    dPv_dot, ...
    'ACE',        ACE, ...
    'phi_db',     phi_db, ...
    'db_active',  db_active, ...
    'grc_active', grc_active, ...
    'lim_active', [sig_Pv < 1, sig_Pm < 1], ...
    'aw_sigma',   [sig_Pv, sig_Pm], ...
    'COP_eff',    COP_eff, ...
    'Tin_abs',    Tin_abs, ...
    'Tout_abs',   Tout_abs );

end % cpmg_dynamics

%% ==================================================================== %%
%  Local nonlinearity primitives
%% ==================================================================== %%

function y = deadbandHard(u, w)
%DEADBANDHARD  Physical governor dead zone of manuscript Eq.(2).
%   Piecewise linear, continuous, non-differentiable at +/- w.
if u > w
    y = u - w;
elseif u < -w
    y = u + w;
else
    y = 0;
end
end

function y = deadbandSmooth(u, w, k)
%DEADBANDSMOOTH  C^inf surrogate of the dead zone, used only for analysis.
%   Recovers deadbandHard as k -> Inf.  Provided so that the Lipschitz-
%   constant estimation in estimate_lipschitz.m can operate on a smooth field.
y = u - w*tanh(k*u/max(w,eps));
end

function Phi = socSaturation(S, S_lo, S_hi, k)
%SOCSATURATION  Two-sided smooth SoC derating of the BESS control channel.
%
%   Phi(S) = 0.25*(1 + tanh(k*(S - S_lo)))*(1 - tanh(k*(S - S_hi)))
%
%   Phi -> 1 for S_lo << S << S_hi and Phi -> 0 at BOTH bounds.  The legacy
%   single-sided form 0.5*(1 + tanh(k*(S - S_mid))) derated only near the
%   lower bound and permitted unbounded charging at full SoC.
%
%   Phi multiplies u_B in G(x), so Phi -> 0 collapses the BESS column of
%   LgLf_h and the safe set becomes uncontrollable through the fast actuator.
%   This is the mechanism analysed in Proposition 1.
Phi = 0.25*(1 + tanh(k*(S - S_lo)))*(1 - tanh(k*(S - S_hi)));
Phi = min(max(Phi, 0), 1);
end

function [xd, sigma] = clampRate(xv, xd, lo, hi, band)
%CLAMPRATE  Lipschitz anti-windup clamp.
%
%   Scales an outward-pointing derivative to zero over a narrow band of width
%   `band` adjacent to each magnitude limit, instead of switching it off at
%   the limit.  Inward-pointing derivatives are never scaled.
%
%   WHY THIS IS NOT A HARD SWITCH
%   -----------------------------
%   The hard version used previously made F DISCONTINUOUS: at dPv = Pv_max
%   with a positive raw rate, F(4) jumped from +2.425 to 0.  Assumption 1
%   asserts that F is Lipschitz on the compact domain D, and a jump makes
%   that assertion false.  The defect was detected by estimate_lipschitz,
%   which returned a secant estimate LARGER than the supremum of the local
%   Jacobian norm (207.72 against 127.81); a secant ratio can only exceed the
%   Jacobian supremum where the field fails to be differentiable, and can only
%   grow without bound where it fails to be continuous.
%
%   With the banded form the field is continuous and piecewise linear, so it
%   IS Lipschitz, with a local slope of order |xd|/band.  For band = 0.02 and
%   a raw rate of 2.425 p.u./s the induced slope is about 121, comparable to
%   the governor mode 1/(R*Tg) = 125 already present, so Assumption 1 becomes
%   literally true without inflating the constants.
%
%   The clamp still returns EXACTLY zero at the limit itself, so the physical
%   bound is enforced exactly and the verification tests are unaffected.
%
%   sigma in [0,1] is the applied scaling, returned for diagnostics.
if nargin < 5 || isempty(band), band = 1e-3; end
band = max(band, eps);

if xd > 0
    sigma = satlin((hi - xv)/band);
elseif xd < 0
    sigma = satlin((xv - lo)/band);
else
    sigma = 1;
end
xd = sigma*xd;
end

function z = satlin(z)
%SATLIN  Unit ramp saturated to [0,1].
z = min(max(z,0),1);
end
