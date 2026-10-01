function [xn, info] = cpmg_integrate(x, ug, u, d, p, tie_in, h, method)
%CPMG_INTEGRATE  One fixed-step integration of the single-area CPMG plant.
%
%   [xn, info] = CPMG_INTEGRATE(x, ug, u, d, p, tie_in, h, method) advances the
%   state by one step of size h using the requested method.  The edge command
%   u, the AGC command ug, the disturbance d and the tie-line coupling tie_in
%   are held constant across the step (zero-order hold), which is the correct
%   representation of a sampled-data edge controller.
%
%   INPUTS
%     x       9x1   state at t
%     ug      1x1   held AGC command                                [p.u.]
%     u       2x1   held edge command [u_B ; u_A]                   [p.u.]
%     d       3x1   held disturbance                                [-]
%     p       1x1   parameter struct
%     tie_in  1x1   held tie-line coupling term (0 if isolated)     [p.u./s]
%     h       1x1   step size                                       [s]
%     method  char  'rk4' (default) | 'euler'
%
%   OUTPUTS
%     xn      9x1   state at t + h
%     info    struct with the diagnostics from the FIRST stage evaluation
%             (the one at the current state), plus .nfev
%
%   WHY ZOH ON tie_in
%   -----------------
%   In the multi-area runner the coupling depends on the frequencies of the
%   other areas.  Holding it across the RK4 stages makes each area's update
%   independent within a step, which is what a genuinely decentralised
%   implementation can do.  The resulting splitting error is O(h^2) and is
%   quantified by verify_plant.m.  Use cpmg_integrate_multiarea for a fully
%   coupled update when that error must be eliminated.
%
%   See also CPMG_DYNAMICS, CPMG_INTEGRATE_MULTIAREA, VERIFY_PLANT.

if nargin < 8 || isempty(method), method = p.spec.num.integrator; end
if nargin < 6 || isempty(tie_in), tie_in = 0; end

f = @(xx) fieldEval(xx, ug, u, d, p, tie_in);

switch lower(method)
    case 'euler'
        [k1, dg] = f(x);
        xn = x + h*k1;
        info = dg;  info.nfev = 1;

    case 'rk4'
        [k1, dg] = f(x);
        k2 = f(x + 0.5*h*k1);
        k3 = f(x + 0.5*h*k2);
        k4 = f(x +     h*k3);
        xn = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        info = dg;  info.nfev = 4;

    otherwise
        error('cpmg_integrate:method','Unknown integrator "%s"', method);
end

end

%% ==================================================================== %%
function [xdot, dg] = fieldEval(x, ug, u, d, p, tie_in)
[F, G, dg] = cpmg_dynamics(x, ug, d, p, tie_in);
xdot = F + G*u;
end
