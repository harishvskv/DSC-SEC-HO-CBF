function [Xn, DG] = cpmg_integrate_multiarea(X, Ug, U, D, P, h, method)
%CPMG_INTEGRATE_MULTIAREA  One fixed-step integration of the stacked N-area
%   CPMG plant with FULLY COUPLED stages.
%
%   [Xn, DG] = CPMG_INTEGRATE_MULTIAREA(X, Ug, U, D, P, h, method) advances the
%   whole interconnected system by one step.  Unlike calling cpmg_integrate
%   per area, the tie-line coupling is re-evaluated at every RK4 stage, so no
%   operator-splitting error is introduced.
%
%   This is the integrator used for all reported multi-area results.  The
%   split version exists only so that verify_plant can quantify how large the
%   splitting error would have been.
%
%   INPUTS
%     X       (N*9)x1     stacked state
%     Ug       N x1       held AGC commands
%     U       (N*2)x1     held edge commands
%     D        N x 3      held disturbances, row per area
%     P        1 x N      parameter struct array (params_39bus)
%     h        1x1        step size                                   [s]
%     method  char        'rk4' (default) | 'euler'
%
%   See also CPMG_DYNAMICS_MULTIAREA, CPMG_INTEGRATE.

if nargin < 7 || isempty(method), method = P(1).spec.num.integrator; end

f = @(XX) stackedField(XX, Ug, U, D, P);

switch lower(method)
    case 'euler'
        [k1, DG] = f(X);
        Xn = X + h*k1;

    case 'rk4'
        [k1, DG] = f(X);
        k2 = f(X + 0.5*h*k1);
        k3 = f(X + 0.5*h*k2);
        k4 = f(X +     h*k3);
        Xn = X + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

    otherwise
        error('cpmg_integrate_multiarea:method','Unknown integrator "%s"', method);
end

end

%% ==================================================================== %%
function [Xdot, DG] = stackedField(X, Ug, U, D, P)
[F, G, DG] = cpmg_dynamics_multiarea(X, Ug, D, P);
Xdot = F + G*U;
end
