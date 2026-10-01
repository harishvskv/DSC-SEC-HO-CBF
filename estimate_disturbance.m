function [d_hat, info] = estimate_disturbance(x, rocof, p)
%ESTIMATE_DISTURBANCE  Reconstruct the unmeasured net power disturbance from
%   measured RoCoF and measured states, by inverting the swing equation.
%
%   [d_hat, info] = ESTIMATE_DISTURBANCE(x, rocof, p)
%
%   WHY THIS IS NEEDED
%   ------------------
%   The barrier derivative Lf h = d(df)/dt depends on the renewable injection
%   and the load, neither of which the edge controller measures.  Evaluating
%   the barrier with a zero disturbance therefore reports an APPARENTLY HEALTHY
%   margin rate while the frequency is in fact collapsing.  In the Layer 5
%   gate this made the predictive trigger fire at exactly the instant the
%   barrier was lost, giving a timing margin of precisely zero.  The same
%   defect would have silently weakened the HO-CBF-QP in Layer 6, which builds
%   its constraint from the same Lie derivatives.
%
%   THIS IS NOT AN OBSERVER
%   -----------------------
%   No unmeasured STATE is reconstructed and no dynamics are integrated.  The
%   swing equation is inverted algebraically for the single unmeasured input,
%   using quantities already established as locally available in the edge
%   measurability table:
%
%       M*df_dot = dPm + dPB + dPres - dPA - dPtie - (P_L(df,dPL) - PL0)
%
%   Every term except the disturbance group is measured, and df_dot is
%   measured directly as RoCoF.  Solving for the group gives
%
%       w = M*rocof - [ dPm + dPB - dPA - dPtie - (P_L(df,0) - PL0) ]
%
%   w is the NET unmeasured injection, lumping renewable output and the load
%   step into one scalar.  Separating them is neither possible nor necessary:
%   only their sum enters the barrier.
%
%   RETURNED IN THE dPres CHANNEL
%   -----------------------------
%   d_hat places w in the dP_res slot and leaves dP_L and dT_out at zero.  By
%   construction, calling cpmg_dynamics or barrier_functions with d_hat
%   reproduces the measured RoCoF exactly, so Lf h is then correct.  The second
%   derivative Lf^2 h assumes w is held constant across the sampling interval,
%   which is the same zero-order-hold assumption already made for the control.
%   The resulting mismatch is bounded by the disturbance rate of change and is
%   covered by the tightening term Delta_cbf of Theorem 2.
%
%   PRACTICAL NOTE
%   --------------
%   RoCoF must be filtered in practice.  The filter lag adds directly to the
%   effective age of data and should be included in the timing budget alongside
%   T_pkt and the WAN latency; p.etm.T_rocof records it.
%
%   INPUTS
%     x      9x1  canonical state, all entries locally measured
%     rocof  1x1  measured d(df)/dt [p.u./s]
%     p      parameter struct
%
%   OUTPUTS
%     d_hat  3x1  [w ; 0 ; 0], suitable for cpmg_dynamics / barrier_functions
%     info   struct .w .measured_terms .PL_local .rocof
%
%   See also DA_ETM, BARRIER_FUNCTIONS, CPMG_DYNAMICS.

spec = p.spec;  ix = spec.ix;
eps0 = spec.num.eps_guard;

df    = x(ix.df);
dPtie = x(ix.dPtie);
dPm   = x(ix.dPm);
dPB   = x(ix.dPB);
dPA   = x(ix.dPA);

% Load evaluated at the MEASURED frequency with no load step, which is the
% part of the composite load the edge can compute from local data alone.
if spec.nl.load_frequency
    base    = max(1 + p.D*df, eps0);
    PL_loc  = p.PL0*base^p.gamma;
else
    PL_loc  = p.PL0 + p.D*df;
end

measured = dPm + dPB - dPA - dPtie - (PL_loc - p.PL0);

w = p.M*rocof - measured;

d_hat = zeros(spec.nd,1);
d_hat(spec.id.dPres) = w;

info = struct('w',w,'measured_terms',measured,'PL_local',PL_loc,'rocof',rocof);
end
