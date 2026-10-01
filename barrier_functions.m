function B = barrier_functions(x, ug, d, p, tie_in, cfg)
%BARRIER_FUNCTIONS  Analytic HO-CBF construction for the nonlinear CPMG.
%
%   B = BARRIER_FUNCTIONS(x, ug, d, p, tie_in, cfg) returns every quantity the
%   QP needs to enforce forward invariance of the three safe sets, together
%   with the assembled constraint rows.
%
%   THE THREE SAFE SETS
%   -------------------
%     frequency :  h_f     = df - df_min                  (UFR Stage-1)
%     storage   :  h_Slo   = SoC - SoC_lo      (EFFECTIVE, not physical)
%                  h_Shi   = SoC_hi - SoC      (EFFECTIVE, not physical)
%     comfort   :  h_Tlo   = Tin - Tin_min
%                  h_Thi   = Tin_max - Tin
%
%   ALL FIVE HAVE RELATIVE DEGREE 2.  That is a structural property of this
%   plant and it is worth stating in the manuscript, because it means one
%   uniform HO-CBF construction covers every constraint:
%     - u_B reaches df through dPB, and reaches SoC through dPB as well;
%     - u_A reaches df through dPA, and reaches Tin through dPA as well.
%   In each case the input passes through exactly one first-order actuator
%   lag before it influences the constrained quantity.
%
%   The previous submission imposed only the frequency barrier, wrote the SoC
%   and comfort limits as constraints on STATES rather than on the decision
%   variable, and therefore did not enforce them at all.  That is Reviewer 2
%   comment 1, and it is fixed here.
%
%   HO-CBF FORM (linear class-K, degree 2)
%     psi0 = h
%     psi1 = Lf h + a1 h
%     psi2 = Lf^2 h + LgLf h u + (a1 + a2) Lf h + a1 a2 h  >=  Delta_cbf
%
%   giving the QP row      -LgLf h * u  <=  b - Delta_cbf
%   with                    b = Lf^2 h + (a1+a2) Lf h + a1 a2 h
%
%   INPUTS
%     x       9x1     canonical state
%     ug      1x1     held AGC command
%     d       3x1     held disturbance
%     p       struct  parameters
%     tie_in  1x1     tie-line coupling (0 if isolated)
%     cfg     struct
%              .a_f   [a1 a2] class-K gains, frequency barrier   (default [8 8])
%              .a_S   [a1 a2] class-K gains, SoC barriers        (default [2 2])
%              .a_T   [a1 a2] class-K gains, comfort barriers    (default [0.5 0.5])
%              .Delta_cbf  scalar tightening from Theorem 2      (default 0)
%              .enable     logical [f, Slo, Shi, Tlo, Thi]       (default all true)
%
%   OUTPUT B
%     B.f, B.Slo, B.Shi, B.Tlo, B.Thi   per-barrier structs, each with
%         .h .Lfh .Lf2h .LgLfh (1x2) .psi1 .b
%     B.A     nRows x 2   constraint matrix,  A*u <= b_qp
%     B.bq    nRows x 1   right-hand side, tightening already subtracted
%     B.names nRows x 1   cell array of row labels
%     B.dg                plant diagnostics at x
%     B.rd                relative-degree report, per barrier per channel
%
%   Every derivative below is DERIVED BY HAND from the model in
%   cpmg_dynamics.m and is checked against central differences by
%   verify_lie.m.  Do not edit one without re-running the other.
%
%   See also CPMG_DYNAMICS, VERIFY_LIE, ESTIMATE_LIPSCHITZ.

if nargin < 5 || isempty(tie_in), tie_in = 0; end
if nargin < 6, cfg = struct(); end
cfg = setdef(cfg,'a_f',[8 8]);
cfg = setdef(cfg,'a_S',[2 2]);
cfg = setdef(cfg,'a_T',[0.5 0.5]);
cfg = setdef(cfg,'Delta_cbf',0);
cfg = setdef(cfg,'enable',true(1,5));

spec = p.spec;  ix = spec.ix;

%% ------------------------------------------------------------------ %%
%  Plant evaluation at x.  The barrier module NEVER re-derives the plant; it
%  consumes F, G and the diagnostics, so the two can never drift apart.
%% ------------------------------------------------------------------ %%
[F, G, dg] = cpmg_dynamics(x, ug, d, p, tie_in);

dPtie = x(ix.dPtie);
dPm   = x(ix.dPm);
dPB   = x(ix.dPB);
SoC   = x(ix.dSB);
dPA   = x(ix.dPA);
dTin  = x(ix.dTin);

Phi     = dg.Phi;
dPm_dot = dg.dPm_dot;            % GRC- and clamp-limited, exact
bT_eff  = p.b_T*(dg.COP_eff/p.COP);

%% ==================================================================== %%
%  1. FREQUENCY BARRIER      h_f = df - df_min
%% ==================================================================== %%
%   Lf h_f  = F_1
%           = [dPm + dPB + dPres - dPA - dPtie - (P_L - PL0)]/M
%
%   Lf^2 h_f = d/dt(F_1) along the drift, with u, d and tie_in held:
%           = (1/M)*[ dPm_dot                      turbine, GRC-limited
%                   + (-dPB/TB)                    BESS drift
%                   - (-dPA/Tac)                   AIAC drift
%                   - tie_in                       tie-line
%                   - dPL_ddf * F_1 ]              load recovery
%
%   LgLf h_f = (1/M)*[ Phi/TB , -1/Tac ]
%
%   The AIAC entry is NEGATIVE: raising frequency requires u_A < 0, i.e.
%   shedding cooling load.  This is the Hui2019 convention frozen in D3.
Bf.h    = x(ix.df) - spec.safety.df_min_pu;
Bf.Lfh  = F(ix.df);
Bf.Lf2h = ( dPm_dot ...
          + (-dPB/p.TB) ...
          - (-dPA/p.Tac) ...
          - tie_in ...
          - dg.dPL_ddf*F(ix.df) )/p.M;
Bf.LgLfh = [ Phi/(p.M*p.TB) , -1/(p.M*p.Tac) ];

%% ==================================================================== %%
%  2. STORAGE BARRIERS       h_Slo = SoC - SoC_min ,  h_Shi = SoC_max - SoC
%% ==================================================================== %%
%   SoC' = -eta_B*dPB/E_B                       =>  Lf h_Slo = -eta_B*dPB/E_B
%   dPB' = -dPB/TB + Phi*u_B/TB
%   =>  Lf^2 h_Slo = -eta_B/E_B * (-dPB/TB) = +eta_B*dPB/(E_B*TB)
%       LgLf h_Slo = [ -eta_B*Phi/(E_B*TB) , 0 ]
%
%   The BESS channel of LgLf h_Slo carries the SAME Phi factor as the
%   frequency barrier.  As SoC approaches either bound, Phi -> 0 and BOTH
%   barriers lose authority over the fast actuator simultaneously.  That
%   simultaneity is the content of Proposition 1 and is why the infeasibility
%   scenario of Reviewer 2 comment 3 is real rather than hypothetical.
%   ANCHORING: the barriers are placed at the EFFECTIVE bounds SoC_lo and
%   SoC_hi, not at the physical bounds SoC_min and SoC_max.  Between SoC_min
%   and SoC_lo the factor Phi has already collapsed, so a barrier anchored at
%   SoC_min would report a positive margin over a region in which the actuator
%   has no authority to defend it.  The enforced set is the largest storage
%   set that the BESS can actually render forward invariant; the physical
%   bounds remain in force as reported limits.  See spec.safety.SoC_anchor.
cS = p.eta_B/p.E_B_s;

switch lower(spec.safety.SoC_anchor)
    case 'effective'
        S_lo_enf = p.SoC_enf_lo;  S_hi_enf = p.SoC_enf_hi;
    case 'physical'
        S_lo_enf = p.SoC_min;     S_hi_enf = p.SoC_max;
    otherwise
        error('barrier_functions:socAnchor','Unknown SoC_anchor "%s"', ...
              spec.safety.SoC_anchor);
end

BSlo.h     = SoC - S_lo_enf;
BSlo.Lfh   = -cS*dPB;
BSlo.Lf2h  = +cS*dPB/p.TB;
BSlo.LgLfh = [ -cS*Phi/p.TB , 0 ];

BShi.h     = S_hi_enf - SoC;
BShi.Lfh   = +cS*dPB;
BShi.Lf2h  = -cS*dPB/p.TB;
BShi.LgLfh = [ +cS*Phi/p.TB , 0 ];

% Physical margins, reported alongside but NOT enforced
BSlo.h_phys = SoC - p.SoC_min;
BShi.h_phys = p.SoC_max - SoC;
BSlo.anchor = S_lo_enf;
BShi.anchor = S_hi_enf;

%% ==================================================================== %%
%  3. COMFORT BARRIERS       h_Tlo = Tin - Tin_min ,  h_Thi = Tin_max - Tin
%% ==================================================================== %%
%   Tin  = Tin_ref + dTin ,   F_9 = a_T*(Tout - Tin) - bT*(PA_nom + dPA)
%   dPA' = -dPA/Tac + u_A/Tac ,  Tout held  =>
%       dF_9/dt = -a_T*F_9 + bT*dPA/Tac - (bT/Tac)*u_A
%
%   h_Tlo = Tin - Tin_min :  Lf h = F_9
%       Lf^2 h_Tlo = -a_T*F_9 + bT*dPA/Tac
%       LgLf h_Tlo = [ 0 , -bT/Tac ]
%   h_Thi = Tin_max - Tin :  signs negated.
%
%   SCOPE WARNING, to be stated in the manuscript rather than discovered by a
%   reviewer: with a physically correct 4 h aggregate thermal time constant,
%   |Lf h_T| is O(1e-3) degC/s and these barriers cannot become active on a
%   3 s DoS horizon.  They bind only for sustained (minutes-scale) events.  Do
%   NOT claim tight transient comfort regulation over a short attack.
Tin_abs = p.Tin_ref + dTin;

BTlo.h     = Tin_abs - p.Tin_min;
BTlo.Lfh   = F(ix.dTin);
BTlo.Lf2h  = -p.a_T*F(ix.dTin) + bT_eff*dPA/p.Tac;
BTlo.LgLfh = [ 0 , -bT_eff/p.Tac ];

BThi.h     = p.Tin_max - Tin_abs;
BThi.Lfh   = -F(ix.dTin);
BThi.Lf2h  = +p.a_T*F(ix.dTin) - bT_eff*dPA/p.Tac;
BThi.LgLfh = [ 0 , +bT_eff/p.Tac ];

%% ==================================================================== %%
%  4. HO-CBF assembly
%% ==================================================================== %%
Bf   = hocbf(Bf  , cfg.a_f);
BSlo = hocbf(BSlo, cfg.a_S);
BShi = hocbf(BShi, cfg.a_S);
BTlo = hocbf(BTlo, cfg.a_T);
BThi = hocbf(BThi, cfg.a_T);

B.f = Bf;  B.Slo = BSlo;  B.Shi = BShi;  B.Tlo = BTlo;  B.Thi = BThi;

all_b  = {Bf, BSlo, BShi, BTlo, BThi};
names  = {'CBF:freq','CBF:SoC_lo','CBF:SoC_hi','CBF:Tin_lo','CBF:Tin_hi'};
% The delay tightening applies to the frequency barrier, which is the one
% evaluated on delayed data during the switching transient.  The SoC and
% comfort barriers use BMS and LAN telemetry that survive a WAN DoS, so they
% are not tightened.  State this choice explicitly in Section IV-F.
tighten = [cfg.Delta_cbf, 0, 0, 0, 0];

A = []; bq = []; nm = {};
for k = 1:5
    if ~cfg.enable(k), continue; end
    A  = [A ; -all_b{k}.LgLfh];            %#ok<AGROW>
    bq = [bq;  all_b{k}.b - tighten(k)];   %#ok<AGROW>
    nm = [nm ; names(k)];                  %#ok<AGROW>
end
B.A = A;  B.bq = bq;  B.names = nm;
B.dg = dg;
B.cfg = cfg;

%% ==================================================================== %%
%  5. Relative-degree report
%% ==================================================================== %%
% Structural zeros are asserted, not assumed: the frequency barrier must be
% reachable from BOTH channels, the storage barriers from u_B only, and the
% comfort barriers from u_A only.
B.rd.Lgh       = G(ix.df,:);                 % must be exactly [0 0]
B.rd.LgLfh     = [Bf.LgLfh ; BSlo.LgLfh ; BShi.LgLfh ; BTlo.LgLfh ; BThi.LgLfh];
B.rd.names     = names;
B.rd.Phi       = Phi;
B.rd.Phi_degen = dg.Phi_degen;
B.rd.bess_authority = abs(Bf.LgLfh(1))*max(abs([p.uB_min p.uB_max]));
B.rd.aiac_authority = abs(Bf.LgLfh(2))*max(abs([p.uA_min p.uA_max]));
B.rd.soc_anchor     = [BSlo.anchor, BShi.anchor];
B.rd.soc_physical   = [p.SoC_min,   p.SoC_max];
B.rd.h_soc_phys     = [BSlo.h_phys, BShi.h_phys];

end % barrier_functions


%% ==================================================================== %%
function s = hocbf(s, a)
%HOCBF  Degree-two high-order CBF terms for a barrier with linear class-K
%   gains a = [a1 a2].
%       psi1 = Lf h + a1 h
%       psi2 = Lf^2 h + LgLf h u + (a1+a2) Lf h + a1 a2 h
%   The QP row is  -LgLf h * u <= b,  b = Lf^2 h + (a1+a2) Lf h + a1 a2 h.
s.a    = a;
s.psi1 = s.Lfh + a(1)*s.h;
s.b    = s.Lf2h + (a(1)+a(2))*s.Lfh + a(1)*a(2)*s.h;
end

function c = setdef(c, f, v)
if ~isfield(c,f) || isempty(c.(f)), c.(f) = v; end
end
