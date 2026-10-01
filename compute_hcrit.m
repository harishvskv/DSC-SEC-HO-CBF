function [h_crit, info] = compute_hcrit(p, dP, opt)
%COMPUTE_HCRIT  Critical initial margin below which a step imbalance cannot
%   be arrested by ANY controller with the given actuators.
%
%   [h_crit, info] = COMPUTE_HCRIT(p, dP, opt)
%
%   WHAT THIS QUANTIFIES
%   --------------------
%   Forward invariance of the frequency safe set is not a property of the
%   control law alone.  Even a controller that commands full output at the
%   first sample cannot prevent a breach if the actuators need longer to
%   deliver that output than the margin survives.  For a step imbalance dP
%   there is therefore a critical initial margin h_crit such that
%
%       h0 >= h_crit   ->  the barrier can be held
%       h0 <  h_crit   ->  a breach is unavoidable, whatever the controller
%
%   h_crit is fixed by the actuator dynamics: the BESS power-loop time
%   constant T_B, the command slew limit u_B_rate, the authority u_B_max and
%   the inertia M.  It is NOT reduced by better tuning, and reporting it is
%   the honest alternative to adjusting parameters until a test passes.
%
%   ANALYTIC DECOMPOSITION
%   ----------------------
%   With the command driven to u_B_max and the power loop a first-order lag,
%   dP_B(t) = u_B_max*(1 - exp(-t/T_B)), the deficit vanishes at
%
%       t* = -T_B*ln(1 - dP/u_B_max)
%
%   and the margin consumed up to that instant is
%
%       loss_lag = (1/M)*[ (dP - u_B_max)*t* + dP*T_B ]
%
%   The command itself cannot step, so the slew adds approximately
%
%       loss_slew = (dP/(2*M))*(dP/u_B_rate)
%
%   and zero-order-hold sampling adds at most (dP/M)*Ts_ctrl.  The sum is an
%   estimate; the value returned as h_crit is MEASURED by bisection on the
%   closed loop, and the estimate is returned alongside for interpretation.
%
%   THE AIAC DOES NOT APPEAR
%   ------------------------
%   At u_A_rate = 0.10 p.u./s the air-conditioning fleet moves 0.001 p.u. per
%   controller period and needs seconds to reach any useful contribution,
%   while the margin survives tenths of a second.  Its contribution to h_crit
%   is negligible and it is omitted from the estimate rather than included
%   with a coefficient that would imply it helps.
%
%   INPUTS
%     p    parameter struct
%     dP   step imbalance [p.u.]
%     opt  .tol       bisection tolerance on h0 (default 1e-4)
%          .Tend      simulation horizon [s] (default 6)
%          .qp        options forwarded to ctrl_hocbf_qp
%          .verbose   print a summary (default false)
%
%   OUTPUTS
%     h_crit  measured critical margin [p.u.].  Inf when dP exceeds the total
%             actuator authority, in which case no margin is sufficient.
%     info    .est_lag .est_slew .est_zoh .est_total
%             .t_star  time for the BESS to cover the deficit
%             .frac    h_crit as a fraction of the full safe set
%             .authority_exceeded
%
%   See also CTRL_HOCBF_QP, COMPUTE_SCATD, VERIFY_QP.

if nargin < 3, opt = struct(); end
opt = setdef(opt,'tol',1e-4);
opt = setdef(opt,'Tend',6.0);
opt = setdef(opt,'verbose',false);
opt = setdef(opt,'qp',struct());

spec = p.spec;  ix = spec.ix;
Ts   = spec.num.Ts_plant;
Tc   = spec.num.Ts_ctrl;
h_max = -spec.safety.df_min_pu;

%% ------------------------------------------------------------------ %%
%  Analytic estimate
%% ------------------------------------------------------------------ %%
info = struct();
info.authority_exceeded = dP >= p.uB_max;

if info.authority_exceeded
    info.t_star = Inf;
    info.est_lag = Inf; info.est_slew = Inf; info.est_zoh = Inf;
    info.est_total = Inf;
    h_crit = Inf;
    info.frac = Inf;
    if opt.verbose, printH(h_crit, info, p, dP, h_max); end
    return
end

t_star        = -p.TB*log(1 - dP/p.uB_max);
info.t_star   = t_star;
info.est_lag  = ( (dP - p.uB_max)*t_star + dP*p.TB )/p.M;
info.est_slew = (dP/(2*p.M))*(dP/p.uB_rate);
info.est_zoh  = (dP/p.M)*Tc;
info.est_total = info.est_lag + info.est_slew + info.est_zoh;

%% ------------------------------------------------------------------ %%
%  Measured, by bisection on the initial margin
%% ------------------------------------------------------------------ %%
holds = @(h0) barrierHeld(h0, dP, p, opt, Ts, Tc, ix, spec);

lo = 0;  hi = h_max;
if ~holds(hi)
    h_crit = Inf;                       % even a full margin is insufficient
else
    while (hi - lo) > opt.tol
        mid = 0.5*(lo + hi);
        if holds(mid), hi = mid; else, lo = mid; end
    end
    h_crit = hi;
end
info.frac = h_crit/h_max;

if opt.verbose, printH(h_crit, info, p, dP, h_max); end
end % compute_hcrit


%% ==================================================================== %%
function ok = barrierHeld(h0, dP, p, opt, Ts, Tc, ix, spec)
x = p.x0;  x(ix.df) = spec.safety.df_min_pu + h0;
d = [0; dP; 0];
cs = [];  u = [0;0];
cfg0 = struct('a_f',[8 8],'Delta_cbf',0);
nStep = round(opt.Tend/Ts);  kCtrl = round(Tc/Ts);
ok = true;
for n = 1:nStep
    if mod(n-1,kCtrl) == 0
        [~, uref] = ctrl_pi_agc(x, p, 'edge_only');
        F  = cpmg_dynamics(x, 0, d, p, 0);
        [u, cs] = ctrl_hocbf_qp(x, F(ix.df), uref, p, cs, opt.qp);
    end
    B = barrier_functions(x, 0, d, p, 0, cfg0);
    if B.f.h < 0, ok = false; return; end
    x = cpmg_integrate(x, 0, u, d, p, 0, Ts, 'rk4');
end
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

function printH(h_crit, info, p, dP, h_max)
fprintf('\n--- critical initial margin, dP = %.4f p.u. ---\n', dP);
if info.authority_exceeded
    fprintf('  imbalance exceeds u_B_max = %.2f p.u.: no initial margin suffices\n\n', p.uB_max);
    return
end
fprintf('  t* (BESS covers deficit) = %.4f s   [T_B = %.3f s]\n', info.t_star, p.TB);
fprintf('  loss  lag   %.5f\n', info.est_lag);
fprintf('  loss  slew  %.5f   [u_B_rate = %.1f p.u./s]\n', info.est_slew, p.uB_rate);
fprintf('  loss  ZOH   %.5f\n', info.est_zoh);
fprintf('  estimate    %.5f\n', info.est_total);
fprintf('  MEASURED h_crit = %.5f p.u.  (%.1f%% of the safe set %.4f)\n\n', ...
        h_crit, 100*info.frac, h_max);
end
