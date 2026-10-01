function V = validate_scatd(p, L, R, opt)
%VALIDATE_SCATD  Empirical test of the SC-ATD bound.
%
%   V = VALIDATE_SCATD(p, L, R, opt) finds, by bisection on the actual attack
%   duration in closed-loop simulation, the longest DoS the plant survives from
%   each of a grid of initial operating states, and compares that measurement
%   against the predictions of COMPUTE_SCATD.
%
%   WHY THIS IS THE ONLY THING THAT MAKES THEOREM 1 A THEOREM AND NOT A CLAIM
%   ------------------------------------------------------------------------
%   The submitted manuscript reported a single SC-ATD value with no derivation
%   path and no measurement, and then ran a stochastic study whose attack
%   durations all lay INSIDE that value while reporting a 21.2 percent failure
%   rate.  A bound that is violated by a fifth of the samples that satisfy it
%   is not a bound.  This routine makes the prediction falsifiable: it sweeps
%   the initial margin, measures the true tolerable duration, and reports the
%   signed error.  A bound that is ever OPTIMISTIC (predicted longer than
%   measured) is a failure and is flagged as such.
%
%   WHAT IS VALIDATED HERE AND WHAT IS NOT
%   --------------------------------------
%   The UNMITIGATED bound T_quad is fully validated: it needs only the frozen
%   edge command, which requires no controller beyond ctrl_pi_agc.
%
%   The PROPOSED bound T_prop depends on the local HO-CBF-QP, which is Layer 6
%   and does not exist yet.  Its energy-limited component T_sustain is checked
%   here against an OPEN-LOOP surrogate in which the BESS is commanded to
%   supply the imbalance directly.  Full closed-loop validation is deferred to
%   Layer 7 and this routine says so in its report rather than implying the
%   result is complete.
%
%   INPUTS
%     p, L, R  as in compute_scatd
%     opt      optional:
%               .h0_grid    initial margins to sweep       (default 5 values)
%               .dP_imb     imbalance for the sweep        (default 0.10)
%               .SoC_grid   initial SoC values             (default 3 values)
%               .tol        bisection tolerance [s]        (default 1e-3)
%               .Tmax       upper bracket [s]              (default 60)
%               .a_f        class-K gains                  (default [8 8])
%               .verbose    print a table                  (default true)
%
%   OUTPUTS (struct V)
%     V.unmit    table: h0, hdot0, T_measured, T_quad, T_exp, signed error
%     V.energy   table: SoC0, T_measured, T_sustain, signed error
%     V.worst_optimism   largest amount by which any prediction EXCEEDS the
%                        measurement.  Must be <= 0 for the bound to hold.
%     V.pass     true when no prediction is optimistic
%
%   See also COMPUTE_SCATD, ESTIMATE_CONTRACTION_RATES.

if nargin < 4, opt = struct(); end
opt = setdef(opt,'h0_grid',[0.003 0.005 0.007 0.009 0.011]);
% Default imbalance is set ABOVE the critical value beta*|df_min|, because
% below it the governor arrests the decline and the barrier is never
% approached, so the validation would be vacuous.
opt = setdef(opt,'dP_imb', 1.2*p.beta*abs(p.spec.safety.df_min_pu));
opt = setdef(opt,'SoC_grid',[0.35 0.50 0.70]);
opt = setdef(opt,'tol',1e-3);
opt = setdef(opt,'Tmax',60);
opt = setdef(opt,'a_f',[8 8]);
opt = setdef(opt,'verbose',true);

spec = p.spec;  ix = spec.ix;
cfg  = struct('a_f',opt.a_f,'Delta_cbf',0);

%% ==================================================================== %%
%  A.  UNMITIGATED BOUND.  Frozen edge command, AGC blocked.
%% ==================================================================== %%
U = struct('h0',{},'hdot0',{},'h_min',{},'T_meas',{},'T_quad',{}, ...
           'T_quad_fit',{},'T_exp',{},'err_quad',{},'err_fit',{},'optimistic',{});

for k = 1:numel(opt.h0_grid)
    h0 = opt.h0_grid(k);

    x0 = p.x0;
    x0(ix.df) = spec.safety.df_min_pu + h0;

    st = struct('x0',x0,'SoC0',p.SoC_nom,'dP_imb',opt.dP_imb);
    S  = compute_scatd(p, L, R, st, struct('a_f',opt.a_f,'use_bound',true));

    % Measured: longest attack for which h stays at or above the floor
    survive = @(T) runFrozen(x0, T, opt.dP_imb, p, cfg, S.h_min);
    T_meas  = bisectSurvival(survive, 0, opt.Tmax, opt.tol);

    % A prediction is OPTIMISTIC when it exceeds the measurement.  Only the
    % rigorous bound is allowed to fail the gate; the fitted value is
    % reported for accuracy but never supports a safety claim.
    err  = S.T_quad_bound - T_meas;
    errf = S.T_quad_fit   - T_meas;
    U(end+1) = struct('h0',h0,'hdot0',S.state.hdot0,'h_min',S.h_min, ...
                      'T_meas',T_meas,'T_quad',S.T_quad_bound, ...
                      'T_quad_fit',S.T_quad_fit,'T_exp',S.T_exp, ...
                      'err_quad',err,'err_fit',errf, ...
                      'optimistic',err > opt.tol); %#ok<AGROW>
end
V.unmit = U;

%% ==================================================================== %%
%  B.  ENERGY-LIMITED SUSTAIN TIME, open-loop surrogate.
%% ==================================================================== %%
%   The BESS is commanded to supply the imbalance directly.  This bounds the
%   closed-loop behaviour from one side: the real QP will command no more
%   than this to hold the barrier, so the measured drain time is a lower
%   bound on what the closed loop achieves.
E = struct('SoC0',{},'T_meas',{},'T_sustain',{},'err',{},'optimistic',{});

for k = 1:numel(opt.SoC_grid)
    SoC0 = opt.SoC_grid(k);
    xs = p.x0;  xs(ix.dSB) = SoC0;
    st = struct('x0',xs,'SoC0',SoC0,'dP_imb',opt.dP_imb);
    S  = compute_scatd(p, L, R, st, struct('a_f',opt.a_f,'use_bound',true));

    x0 = p.x0;  x0(ix.dSB) = SoC0;
    % Adaptive horizon: the analytic drain time plus generous margin, so the
    % measurement is never truncated by the cap.
    Tcap   = 2*max(S.T_sustain, 1) + 60;
    T_meas = drainTime(x0, S.P_bess, p, Tcap);

    err = S.T_sustain - T_meas;
    E(end+1) = struct('SoC0',SoC0,'T_meas',T_meas,'T_sustain',S.T_sustain, ...
                      'err',err,'optimistic',err > 10*opt.tol); %#ok<AGROW>
end
V.energy = E;

%% ==================================================================== %%
V.worst_optimism = max([ [U.err_quad] , [E.err] ]);
V.pass           = ~any([U.optimistic]) && ~any([E.optimistic]);
V.closed_loop_validated = false;   % requires Layer 6, see header
V.meta = struct('dP_imb',opt.dP_imb,'a_f',opt.a_f,'plant',p.name, ...
                'spec',spec.version);

if opt.verbose, printV(V,p); end
if ~V.pass
    warning('validate_scatd:OPTIMISTIC', ...
        ['The SC-ATD bound is OPTIMISTIC for at least one state. A bound that ' ...
         'predicts longer survival than the plant delivers is not a bound. ' ...
         'Increase the safety inflation in estimate_lipschitz, or use the ' ...
         'analytic a_d_bound rather than the fitted value, before reporting.']);
end
end % validate_scatd


%% ==================================================================== %%
%  Simulation kernels
%% ==================================================================== %%

function ok = runFrozen(x0, T_attack, dP, p, cfg, h_floor)
%RUNFROZEN  True when the barrier stays at or above h_floor for the whole
%   attack, with the edge command frozen and the AGC blocked.
spec  = p.spec;
hstep = spec.num.Ts_plant;
d     = [0; dP; 0];

[~, ue] = ctrl_pi_agc(x0, p, 'agc_edge');     % last command before the attack
x  = x0;
N  = max(round(T_attack/hstep), 1);
ok = true;
for n = 1:N
    B = barrier_functions(x, 0, d, p, 0, cfg);
    if B.f.h < h_floor, ok = false; return; end
    x = cpmg_integrate(x, 0, ue, d, p, 0, hstep, 'rk4');
end
B = barrier_functions(x, 0, d, p, 0, cfg);
ok = B.f.h >= h_floor;
end

function T = drainTime(x0, P_bess, p, Tmax)
%DRAINTIME  Time for the SoC to fall from its initial value to the ENFORCED
%   lower bound while the BESS supplies P_bess.
spec  = p.spec;  ix = spec.ix;
hstep = spec.num.Ts_plant;
d     = zeros(spec.nd,1);
x = x0;  T = Inf;
N = round(Tmax/hstep);
for n = 1:N
    if x(ix.dSB) <= p.SoC_enf_lo
        T = (n-1)*hstep; return
    end
    % Command the BESS to hold P_bess; Phi scales what actually flows.
    uB = P_bess;
    x  = cpmg_integrate(x, 0, [uB;0], d, p, 0, hstep, 'rk4');
end
end

function T = bisectSurvival(survive, lo, hi, tol)
%BISECTSURVIVAL  Largest T for which survive(T) is true, assuming survival is
%   monotone decreasing in T.
if ~survive(lo), T = 0; return; end
if survive(hi),  T = hi; return; end
while (hi - lo) > tol
    mid = 0.5*(lo + hi);
    if survive(mid), lo = mid; else, hi = mid; end
end
T = lo;
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

function printV(V,p)
fprintf('\n');
fprintf('========================================================================\n');
fprintf(' SC-ATD VALIDATION   plant: %s   imbalance %.2f p.u.\n', ...
        p.name, V.meta.dP_imb);
fprintf('========================================================================\n');
fprintf(' A. UNMITIGATED (frozen edge command, no DA-ETM)\n');
fprintf(' %8s %9s %8s %9s %10s %9s %9s %s\n', ...
        'h0','hdot0','h_min','T_meas','T_qd(bnd)','T_qd(fit)','err','verdict');
fprintf('------------------------------------------------------------------------\n');
for k = 1:numel(V.unmit)
    r = V.unmit(k);
    fprintf(' %8.4f %9.4f %8.5f %9.3f %10.3f %9.3f %+9.3f %s\n', ...
        r.h0, r.hdot0, r.h_min, r.T_meas, r.T_quad, r.T_quad_fit, r.err_quad, ...
        tern(~r.optimistic,'ok','OPTIMISTIC'));
end
fprintf('\n The exponential model predicts T_exp = ');
if all(isinf([V.unmit.T_exp]))
    fprintf('Inf for every state: an exponential\n envelope never reaches the floor, so it cannot bound a UFLS crossing.\n');
else
    fprintf('%s\n', mat2str(round([V.unmit.T_exp],3)));
end

fprintf('\n B. ENERGY-LIMITED SUSTAIN TIME (open-loop surrogate)\n');
fprintf(' %8s %12s %12s %10s %s\n','SoC0','T_meas','T_sustain','err','verdict');
fprintf('------------------------------------------------------------------------\n');
for k = 1:numel(V.energy)
    r = V.energy(k);
    fprintf(' %8.2f %12.3f %12.3f %+10.3f %s\n', ...
        r.SoC0, r.T_meas, r.T_sustain, r.err, tern(~r.optimistic,'ok','OPTIMISTIC'));
end
fprintf('------------------------------------------------------------------------\n');
fprintf(' worst optimism (must be <= 0): %+.4f s\n', V.worst_optimism);
fprintf(' OVERALL: %s\n', tern(V.pass,'PASS','FAIL'));
fprintf(' Closed-loop SC-ATD validation: DEFERRED to Layer 7 (needs the QP).\n');
fprintf('========================================================================\n\n');
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
