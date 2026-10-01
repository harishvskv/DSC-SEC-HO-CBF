function S = compute_scatd(p, L, R, st, opt)
%COMPUTE_SCATD  Safety-Critical Allowable Time Delay, Theorem 1.
%
%   S = COMPUTE_SCATD(p, L, R, st, opt) evaluates the tolerable DoS duration
%   for a given operating state, under three models, and returns the
%   feasibility floor on which all of them rest.
%
%   THE STRUCTURE OF THE BOUND
%   --------------------------
%   Under the corrected switching law the edge controller runs on FRESH LOCAL
%   measurements once the DA-ETM has isolated the WAN.  Stale data therefore
%   matters only during the detection window, which is bounded by tau_max from
%   Proposition 1.  What eventually ends safety is not delay but ENERGY: the
%   BESS drains, Phi collapses, the barrier loses its fast actuator and the QP
%   becomes infeasible.  Hence
%
%       SC-ATD_proposed  =  T_detect  +  T_sustain(h0, SoC0, dP_imb)
%
%   The unmitigated case, in which the edge command stays frozen because there
%   is no DA-ETM, is a separate and much shorter bound and is the one to which
%   the finite-time model applies.
%
%   THREE MODELS, REPORTED SIDE BY SIDE
%     T_exp   legacy exponential envelope, h(t) >= h0*exp(-lambda_d*t).
%             Strictly positive for all finite t, so it cannot predict a
%             crossing.  Reported to show why it is the wrong ansatz.
%     T_quad  finite-time model, h(t) >= h0 + hdot0*t - 0.5*a_d*t^2.
%             This is the unmitigated SC-ATD.
%     T_prop  T_detect + T_sustain, the proposed framework's SC-ATD.
%
%   INPUTS
%     p    parameter struct
%     L    output of estimate_lipschitz
%     R    output of estimate_contraction_rates
%     st   operating state, struct.  EITHER supply a full state
%           .x0      9x1 canonical state at DoS onset  (preferred; every
%                    barrier quantity is then evaluated locally from it)
%           OR the scalars
%           .h0      frequency margin at DoS onset            [p.u.]
%           .hdot0   margin rate at DoS onset, negative       [p.u./s]
%                    (omit to derive it from .dP_imb via R.hdot0)
%           .SoC0    state of charge at DoS onset             [fraction]
%           .dP_imb  standing power imbalance during the attack [p.u.]
%     opt  optional: .a_f class-K gains (default [8 8]);
%                    .use_bound true to use R.a_d_bound instead of R.a_d_fit;
%                    .aiac_share fraction of the imbalance carried by the AIAC;
%                    .tau_local  delay seen by the LOCAL loop after isolation
%                                (default spec.num.Ts_ctrl, NOT tau_max)
%
%   OUTPUTS (struct S)
%     S.h_min       feasibility floor from Proposition 1        [p.u.]
%     S.T_detect    detection window, = L.tau_max_feasible      [s]
%     S.T_exp       legacy exponential bound                    [s]
%     S.T_quad      unmitigated finite-time bound               [s]
%     S.T_sustain   energy-limited sustain time                 [s]
%     S.T_prop      proposed SC-ATD = T_detect + T_sustain      [s]
%     S.binding     which mechanism binds: 'energy' | 'authority' | 'margin'
%     S.P_bess      sustained BESS power assumed                [p.u.]
%     S.detail      intermediate quantities for the appendix
%
%   See also ESTIMATE_CONTRACTION_RATES, VALIDATE_SCATD, ESTIMATE_LIPSCHITZ.

if nargin < 5, opt = struct(); end
opt = setdef(opt,'a_f',[8 8]);
opt = setdef(opt,'use_bound',true);
opt = setdef(opt,'aiac_share',0.0);

spec = p.spec;  ix = spec.ix;
a1 = opt.a_f(1);  a2 = opt.a_f(2);
opt = setdef(opt,'tau_local', spec.num.Ts_ctrl);

st = setdef(st,'dP_imb', p.beta*abs(spec.safety.df_min_pu));

if isfield(st,'x0') && ~isempty(st.x0)
    % Preferred path: evaluate every barrier quantity LOCALLY at the state.
    cfgb = struct('a_f',opt.a_f,'Delta_cbf',0);
    Bx   = barrier_functions(st.x0, 0, [0; st.dP_imb; 0], p, 0, cfgb);
    st.h0     = Bx.f.h;
    st.hdot0  = Bx.f.Lfh;
    st.Lf2h0  = Bx.f.Lf2h;
    st.SoC0   = st.x0(ix.dSB);
else
    st = setdef(st,'SoC0', p.SoC_nom);
    st = setdef(st,'h0',  -spec.safety.df_min_pu);
    if ~isfield(st,'hdot0') || isempty(st.hdot0)
        st.hdot0 = R.hdot0(st.dP_imb);
    end
    st = setdef(st,'Lf2h0', 0);
end

% Two curvature values, reported side by side.
%   a_d_bound is the supremum of -d2h/dt2 over the whole operating domain: a
%     genuine bound, and conservative.
%   a_d_fit   is fitted to the actual DoS trajectories: accurate, but
%     empirical, and it must never be used to support a safety CLAIM.
% Larger a_d always gives a shorter predicted duration, for either sign, so
% max() is the conservative selection in both cases.
a_d_bound = R.a_d_bound;
a_d_fit   = R.a_d_fit;
a_d = a_d_bound;
if ~opt.use_bound, a_d = a_d_fit; end

%% ==================================================================== %%
%  1. FEASIBILITY FLOOR  h_min   (Proposition 1)
%% ==================================================================== %%
%   The tightened HO-CBF row is  -LgLf h * u <= b - Delta_cbf  with u in a box.
%   It admits a solution iff the best achievable barrier authority covers the
%   deficit:
%
%       A_auth + b - Delta_cbf  >=  0,
%       A_auth = sum_j |LgLf h_j| * u_max_j,
%       b      = Lf^2 h + (a1+a2)*Lf h + a1*a2*h
%
%   Solving for h gives the smallest margin at which the QP is still solvable.
%   Below h_min the barrier cannot be enforced no matter what the controller
%   does, so h_min is the correct terminal value for the SC-ATD, NOT zero.
%
%   Note that A_auth carries Phi(SoC0): a depleted battery raises h_min.  The
%   feasibility floor is therefore itself state dependent, which is what makes
%   the whole bound state dependent.
%
%   TWO CHOICES HERE MATTER, and getting either wrong destroys the bound.
%
%   (i)  Lf^2 h and Lf h are taken LOCALLY at the operating state, not at
%        their worst values over the whole operating domain.  Using the
%        domain suprema makes the floor state-INDEPENDENT and, for this
%        plant, larger than the entire safe set, so every predicted
%        tolerable duration collapses to zero.
%
%   (ii) Delta_cbf is evaluated at the delay the LOCAL loop actually sees
%        after isolation, which is the controller period Ts_ctrl, not at
%        tau_max.  tau_max is by construction the delay at which the
%        tightened constraint becomes marginally infeasible, so evaluating
%        the floor there asserts infeasibility by definition.
Phi0 = phiOf(st.SoC0, p);
LgLf = [ Phi0/(p.M*p.TB) , -1/(p.M*p.Tac) ];
A_auth = abs(LgLf(1))*max(abs([p.uB_min p.uB_max])) ...
       + abs(LgLf(2))*max(abs([p.uA_min p.uA_max]));

Delta = L.Delta_cbf_direct(opt.tau_local, a1, a2);

S.h_min = max( ( Delta - A_auth - st.Lf2h0 - (a1+a2)*st.hdot0 ) / (a1*a2), 0 );

% The floor as an explicit function of state of charge.  This curve is the
% quantitative content of Remark 2: a depleted battery has a strictly smaller
% enforceable safe set.  It is the Section VI-E figure.
S.h_min_of_SoC = @(Sx) max( ( Delta ...
    - ( abs(phiOf(Sx,p)/(p.M*p.TB))*max(abs([p.uB_min p.uB_max])) ...
      + abs(1/(p.M*p.Tac))*max(abs([p.uA_min p.uA_max])) ) ...
    - st.Lf2h0 - (a1+a2)*st.hdot0 ) / (a1*a2), 0 );

%% ==================================================================== %%
%  2. DETECTION WINDOW
%% ==================================================================== %%
S.T_detect = L.tau_max_feasible;

%% ==================================================================== %%
%  3. LEGACY EXPONENTIAL BOUND
%% ==================================================================== %%
if isfinite(R.lambda_d) && R.lambda_d > 0 && st.h0 > S.h_min && S.h_min > 0
    S.T_exp = (1/R.lambda_d)*log(st.h0/S.h_min);
else
    S.T_exp = Inf;      % the envelope never reaches the floor
end
S.T_exp_reaches_zero = false;   % an exponential never reaches zero, ever

%% ==================================================================== %%
%  4. UNMITIGATED FINITE-TIME BOUND
%% ==================================================================== %%
%   h(t) = h0 + hdot0*t - 0.5*a_d*t^2 = h_min
%   =>  T = [ hdot0 + sqrt(hdot0^2 + 2*a_d*(h0 - h_min)) ] / a_d
S.T_quad       = time_to_floor(st.h0, st.hdot0, a_d,       S.h_min);
S.T_quad_bound = time_to_floor(st.h0, st.hdot0, a_d_bound, S.h_min);
S.T_quad_fit   = time_to_floor(st.h0, st.hdot0, a_d_fit,   S.h_min);
S.a_d_used     = a_d;

%% ==================================================================== %%
%  5. ENERGY-LIMITED SUSTAIN TIME
%% ==================================================================== %%
%   After isolation the local barrier controller must supply the standing
%   imbalance for as long as the attack lasts.  The AIAC can carry a share of
%   it; the remainder comes from the BESS.  The sustain time is the time until
%   SoC reaches the ENFORCED bound, at which point Phi has degraded to the
%   point where the QP is no longer solvable.
P_aiac = min(opt.aiac_share*st.dP_imb, max(abs([p.uA_min p.uA_max])));
P_bess = max(st.dP_imb - P_aiac, 0);
P_bess = min(P_bess, p.uB_max);
S.P_bess = P_bess;
S.P_aiac = P_aiac;

if P_bess <= 0
    S.T_sustain  = Inf;
    S.binding    = 'authority';
elseif P_bess >= p.uB_max*(1 - 1e-9) && st.dP_imb > p.uB_max + P_aiac
    % The imbalance exceeds the total actuator authority: safety is lost
    % immediately regardless of stored energy.
    S.T_sustain  = 0;
    S.binding    = 'authority';
else
    dSoC        = max(st.SoC0 - p.SoC_enf_lo, 0);
    S.T_sustain = dSoC*p.E_B_s/(p.eta_B*P_bess);
    S.binding   = 'energy';
end

%% ==================================================================== %%
%  6. PROPOSED SC-ATD
%% ==================================================================== %%
S.T_prop = S.T_detect + S.T_sustain;

% If the margin would be lost during the detection window itself, the
% proposed bound collapses to the unmitigated one.
T_margin_only = S.T_quad;
if T_margin_only < S.T_detect
    S.T_prop  = T_margin_only;
    S.binding = 'margin';
end

%% ==================================================================== %%
%  7. Detail
%% ==================================================================== %%
S.detail = struct('Phi0',Phi0,'LgLf',LgLf,'A_auth',A_auth,'Delta_cbf',Delta, ...
                  'a_d',a_d,'a_d_source',ternary(opt.use_bound,'analytic bound','fitted'), ...
                  'h0',st.h0,'hdot0',st.hdot0,'SoC0',st.SoC0,'dP_imb',st.dP_imb, ...
                  'a_f',opt.a_f,'SoC_enf_lo',p.SoC_enf_lo, ...
                  'Lf2h0',st.Lf2h0,'tau_local',opt.tau_local, ...
                  'dP_crit',p.beta*abs(spec.safety.df_min_pu));
S.state = st;
end % compute_scatd


%% ==================================================================== %%
function Phi = phiOf(S, p)
Phi = 0.25*(1 + tanh(p.ks*(S - p.SoC_lo)))*(1 - tanh(p.ks*(S - p.SoC_hi)));
Phi = min(max(Phi,0),1);
end


function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

function o = ternary(c,a,b)
if c, o = a; else, o = b; end
end
