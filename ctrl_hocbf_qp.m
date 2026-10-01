function [u, cs, qi] = ctrl_hocbf_qp(x, rocof, u_ref, p, cs, opt)
%CTRL_HOCBF_QP  Decentralised HO-CBF / CLF quadratic-program safety filter.
%
%   [u, cs, qi] = CTRL_HOCBF_QP(x, rocof, u_ref, p, cs, opt)
%
%   Solves, at each edge controller period,
%
%       min   0.5*(u - u_ref)' W (u - u_ref)  +  0.5*p_delta*delta^2
%      u,delta
%       s.t.  -LgLf h_f    * u  <=  b_f - Delta_cbf        frequency barrier
%             -LgLf h_Slo  * u  <=  b_Slo                  storage, lower
%             -LgLf h_Shi  * u  <=  b_Shi                  storage, upper
%             -LgLf h_Tlo  * u  <=  b_Tlo                  comfort, lower
%             -LgLf h_Thi  * u  <=  b_Thi                  comfort, upper
%              z2*LgLf h_f * u  -  delta  <=  b_clf        CLF, slack relaxed
%              u_min <= u <= u_max,  |u - u_prev| <= rate*Ts
%              delta >= 0
%
%   WHAT THIS FIXES FROM THE SUBMITTED VERSION
%   ------------------------------------------
%   Reviewer 2 comment 1 observed that the framework was presented as CLF based
%   while the optimisation contained no CLF constraint and no slack variable,
%   and that the storage and comfort limits said to be strictly enforced were
%   not imposed on the decision variable at all.  Both are correct.  In the
%   submitted formulation the only row was the frequency barrier; the state of
%   charge appeared nowhere, and the comfort limit was written as a bound on a
%   STATE, which is not an implementable constraint on u.  Here every claimed
%   constraint is a row acting on the decision variable, and the CLF row
%   carries a genuine slack.
%
%   Reviewer 2 comment 3 observed that the SoC-dependent saturation was missing
%   from the barrier constraint.  Phi enters through barrier_functions and is
%   therefore present in both the frequency and the storage rows, with the
%   consequence that both lose fast-actuator authority together.
%
%   THE DISTURBANCE IS RECONSTRUCTED, NOT ASSUMED ZERO
%   -------------------------------------------------
%   Lf h and Lf^2 h depend on the unmeasured net injection.  Evaluating them
%   with a zero disturbance understates the required control action for the
%   whole of a frequency fall while reporting a satisfied constraint, which is
%   the signature of a controller that produces barrier violations from a
%   nominally feasible QP.  rocof is therefore a REQUIRED input and the
%   disturbance is recovered algebraically by estimate_disturbance.
%
%   FEASIBILITY LADDER  (Proposition 1, and the answer to R2-9 and R1-23)
%   --------------------------------------------------------------------
%   The submitted Remark 4 asserted that the solver "inherently executes
%   adaptive constraint relaxation" and that feasibility is "seamlessly
%   restored".  Neither is a property of a QP.  What actually happens is a
%   defined and logged sequence:
%
%     stage 1  full problem
%     stage 2  drop the CLF row: safety outranks performance
%     stage 3  soften the storage and comfort rows with penalised slacks,
%              keeping the frequency barrier hard
%     stage 4  frequency barrier itself infeasible.  Best effort: maximise
%              the barrier derivative subject to actuator limits only, and
%              set qi.safe_guaranteed = false.  NO INVARIANCE CLAIM HOLDS
%              FOR ANY INTERVAL IN WHICH STAGE 4 WAS USED.
%
%   Stage 4 is not a failure of the implementation, it is the physical case
%   in which the imbalance exceeds the available actuator authority.  It must
%   be reported as a fraction of samples, not smoothed over.
%
%   INPUTS
%     x       9x1   local state
%     rocof   1x1   measured d(df)/dt [p.u./s]
%     u_ref   2x1   nominal reference, typically from ctrl_pi_agc 'edge_only'
%     p             parameters
%     cs            controller state, or [] to initialise
%     opt           .Delta_cbf  tightening scalar (default 0)
%                   .a_f .a_S .a_T   class-K gains
%                   .clf        [c1 c2] backstepping gains (default [4 4])
%                   .w          [w_B w_A] cost weights (default [10 500])
%                   .p_delta    CLF slack penalty (default 1e5)
%                   .enable     logical [f Slo Shi Tlo Thi] (default all true)
%                   .use_rate   apply actuator rate rows (default true)
%                   .Ts         controller period (default spec.num.Ts_ctrl)
%
%   OUTPUTS
%     u    2x1   applied edge command
%     cs         updated controller state (.u_prev, counters)
%     qi         diagnostics: .stage .exitflag .delta .psi1 .psi2 .h
%                .active .Phi .sat .safe_guaranteed .solve_time .rows
%
%   The CLF uses a backstepping Lyapunov function on the frequency error,
%       z1 = df,  z2 = df_dot + c1*z1,  V = 0.5*z1^2 + 0.5*z2^2,
%   whose relative degree with respect to u is one, so LgV is nonzero.  A
%   naive V = 0.5*ACE^2 has LgV identically zero here, because u reaches the
%   frequency only at the second derivative, and would give a vacuous row.
%   Area control error regulation is delegated to u_ref, which is computed
%   from locally measured df and dPtie and therefore survives WAN isolation.
%
%   See also BARRIER_FUNCTIONS, ESTIMATE_DISTURBANCE, CTRL_BANK, VERIFY_QP.

spec = p.spec;  ix = spec.ix;

if nargin < 6, opt = struct(); end
opt = setdef(opt,'Delta_cbf',0);
opt = setdef(opt,'a_f',[8 8]);
opt = setdef(opt,'a_S',[2 2]);
opt = setdef(opt,'a_T',[0.5 0.5]);
opt = setdef(opt,'clf',[4 4]);
opt = setdef(opt,'w',[10 500]);
opt = setdef(opt,'p_delta',1e5);
opt = setdef(opt,'enable',true(1,5));
opt = setdef(opt,'use_rate',true);
opt = setdef(opt,'Ts',spec.num.Ts_ctrl);
opt = setdef(opt,'soft_penalty',1e4);

if isempty(cs)
    cs = struct('u_prev',[0;0],'n_stage',zeros(1,4),'n_solve',0);
end

tSolve = tic;

%% ==================================================================== %%
%  1. Barrier evaluation with the reconstructed disturbance
%% ==================================================================== %%
d_hat = estimate_disturbance(x, rocof, p);
cfg   = struct('a_f',opt.a_f,'a_S',opt.a_S,'a_T',opt.a_T, ...
               'Delta_cbf',opt.Delta_cbf,'enable',opt.enable);
B = barrier_functions(x, 0, d_hat, p, 0, cfg);

nRows = size(B.A,1);
Acbf  = [B.A , zeros(nRows,1)];      % barriers do not touch the slack
bcbf  = B.bq;

%% ==================================================================== %%
%  2. CLF row with slack
%% ==================================================================== %%
c1 = opt.clf(1);  c2 = opt.clf(2);
z1 = x(ix.df);
z2 = B.f.Lfh + c1*z1;
V  = 0.5*z1^2 + 0.5*z2^2;

%   Vdot = z1*(z2 - c1*z1) + z2*(Lf2h + LgLfh*u + c1*(z2 - c1*z1))
%   Require Vdot <= -c2*V + delta
Aclf = [ z2*B.f.LgLfh , -1 ];
bclf = -c2*V - z1*(z2 - c1*z1) - z2*( B.f.Lf2h + c1*(z2 - c1*z1) );

%% ==================================================================== %%
%  3. Cost and bounds
%% ==================================================================== %%
W = diag(opt.w);
H = blkdiag(W, opt.p_delta);
f = [ -W*u_ref(:) ; 0 ];
H = (H + H')/2;                       % symmetry guard for quadprog

lb = [ p.uB_min ; p.uA_min ; 0   ];
ub = [ p.uB_max ; p.uA_max ; Inf ];

if opt.use_rate
    % Hui2019 reports an explicit compressor rate limiter.  It is imposed here
    % as a bound on the decision variable rather than inside the plant,
    % because a rate limiter in the plant would destroy control-affineness in
    % u and invalidate the whole HO-CBF construction.
    dmax = [p.uB_rate; p.uA_rate]*opt.Ts;
    lb(1:2) = max(lb(1:2), cs.u_prev - dmax);
    ub(1:2) = min(ub(1:2), cs.u_prev + dmax);
    lb(1:2) = min(lb(1:2), ub(1:2));      % guard against an empty box
end

qopt = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');

%% ==================================================================== %%
%  4. Feasibility ladder
%% ==================================================================== %%
stage = 0;  delta = 0;  soft = zeros(1,4);  exitflag = -99;  z = [];

% ---- stage 1: full problem -------------------------------------------
A1 = [Acbf; Aclf];  b1 = [bcbf; bclf];
[z, ~, exitflag] = quadprog(H, f, A1, b1, [], [], lb, ub, [], qopt);
if exitflag == 1, stage = 1; end

% ---- stage 2: drop the CLF -------------------------------------------
if stage == 0
    [z, ~, exitflag] = quadprog(H, f, Acbf, bcbf, [], [], lb, ub, [], qopt);
    if exitflag == 1, stage = 2; end
end

% ---- stage 3: soften the auxiliary barriers, keep frequency hard ------
if stage == 0
    keep = strcmp(B.names,'CBF:freq');
    nS   = nnz(~keep);
    if nS > 0
        % variables [u; delta; s_1..s_nS]
        Hs = blkdiag(H, opt.soft_penalty*eye(nS));
        fs = [f; zeros(nS,1)];
        Ah = [B.A(keep,:) , zeros(nnz(keep),1) , zeros(nnz(keep),nS)];
        As = [B.A(~keep,:), zeros(nS,1)        , -eye(nS)];
        As2 = [Ah; As];
        bs2 = [B.bq(keep); B.bq(~keep)];
        lbs = [lb; zeros(nS,1)];
        ubs = [ub; inf(nS,1)];
        [zs, ~, exitflag] = quadprog(Hs, fs, As2, bs2, [], [], lbs, ubs, [], qopt);
        if exitflag == 1
            stage = 3;
            z = zs(1:3);
            soft(1:min(nS,4)) = zs(4:min(3+nS,7))';
        end
    end
end

% ---- stage 4: best effort, no invariance claim ------------------------
if stage == 0
    stage = 4;
    % Maximise the barrier derivative LgLf h_f * u over the actuator box.
    % With a linear objective the optimum is at a vertex, obtained directly.
    g = B.f.LgLfh(:);
    ubest = zeros(2,1);
    for j = 1:2
        if g(j) >= 0, ubest(j) = ub(j); else, ubest(j) = lb(j); end
    end
    z = [ubest; 0];
    exitflag = 0;
end

u     = z(1:2);
delta = z(3);

%% ==================================================================== %%
%  5. Diagnostics
%% ==================================================================== %%
psi2 = B.f.b + B.f.LgLfh*u;                       % must be >= Delta_cbf
sat  = [ abs(u(1) - lb(1)) < 1e-9 || abs(u(1) - ub(1)) < 1e-9 , ...
         abs(u(2) - lb(2)) < 1e-9 || abs(u(2) - ub(2)) < 1e-9 ];

qi = struct( ...
    'stage',            stage, ...
    'exitflag',         exitflag, ...
    'delta',            delta, ...
    'soft',             soft, ...
    'h',                B.f.h, ...
    'psi1',             B.f.psi1, ...
    'psi2',             psi2, ...
    'psi2_margin',      psi2 - opt.Delta_cbf, ...
    'V',                V, ...
    'LgLfh',            B.f.LgLfh, ...
    'Phi',              B.dg.Phi, ...
    'Phi_degen',        B.dg.Phi_degen, ...
    'sat',              sat, ...
    'w_hat',            d_hat(spec.id.dPres), ...
    'Delta_cbf',        opt.Delta_cbf, ...
    'safe_guaranteed',  stage <= 3 && psi2 >= opt.Delta_cbf - 1e-9, ...
    'rows',             {B.names}, ...
    'solve_time',       toc(tSolve));

cs.u_prev  = u;
cs.n_stage(stage) = cs.n_stage(stage) + 1;
cs.n_solve = cs.n_solve + 1;
end

%% ==================================================================== %%
function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end
