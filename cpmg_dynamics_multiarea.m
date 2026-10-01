function [F, G, DG] = cpmg_dynamics_multiarea(X, Ug, D, P)
%CPMG_DYNAMICS_MULTIAREA  Stacked nonlinear dynamics of an N-area CPMG.
%
%   [F, G, DG] = CPMG_DYNAMICS_MULTIAREA(X, Ug, D, P) assembles the
%   interconnected plant by calling CPMG_DYNAMICS once per area and adding
%   the tie-line coupling.  Every nonlinearity present in the single-area
%   model (dead-band, GRC, gamma-load, Phi(SoC)) is retained per area; there
%   is no linearised variant of this plant.
%
%   INPUTS
%     X   (N*nx) x 1   stacked state, area i occupying rows (i-1)*nx + (1:nx)
%     Ug   N x 1       central AGC command per area                  [p.u.]
%     D    N x nd      disturbance matrix, row i is area i's [dPres dPL dTout]
%     P    1 x N       struct array from params_39bus, with P(1).net populated
%
%   OUTPUTS
%     F   (N*nx) x 1        stacked drift
%     G   (N*nx) x (N*nu)   block-diagonal edge input matrix.  Block-diagonal
%                           structure is what makes the controller genuinely
%                           DECENTRALISED: area i's inputs cannot influence
%                           area j's states within one integration step except
%                           through the tie-line state.
%     DG  1 x N             per-area diagnostics struct array
%
%   Stacked derivative:   Xdot = F + G*Uedge,  Uedge = [uB_1;uA_1;...;uB_N;uA_N]
%
%   TIE-LINE COUPLING (Assumption 3, declared linear):
%       dPtie_i' = 2*pi * sum_{j ~= i} T_ij * (df_i - df_j)
%
%   See also CPMG_DYNAMICS, PARAMS_39BUS, CPMG_SPEC.

spec = P(1).spec;
nx   = spec.nx;
nu   = spec.nu;
nA   = numel(P);
net  = P(1).net;

assert(numel(X) == nA*nx, 'cpmg_dynamics_multiarea: X must be %d x 1', nA*nx);
assert(numel(Ug) == nA,   'cpmg_dynamics_multiarea: Ug must be %d x 1', nA);
assert(size(D,1) == nA && size(D,2) == spec.nd, ...
       'cpmg_dynamics_multiarea: D must be %d x %d', nA, spec.nd);

F = zeros(nA*nx, 1);
G = zeros(nA*nx, nA*nu);

%% ------------------------------------------------------------------ %%
%  Frequency deviations of every area (needed before the per-area call)
%% ------------------------------------------------------------------ %%
df = zeros(nA,1);
for i = 1:nA
    df(i) = X((i-1)*nx + spec.ix.df);
end

%% ------------------------------------------------------------------ %%
%  Per-area assembly
%% ------------------------------------------------------------------ %%
for i = 1:nA
    rows = (i-1)*nx + (1:nx);
    cols = (i-1)*nu + (1:nu);

    % Tie-line coupling seen by area i
    tie_in = 2*pi * sum( net.T_ij(i,:)' .* (df(i) - df) );

    [Fi, Gi, dgi] = cpmg_dynamics(X(rows), Ug(i), D(i,:)', P(i), tie_in);

    F(rows)       = Fi;
    G(rows, cols) = Gi;

    dgi.tie_in = tie_in;
    dgi.area   = i;
    DG(i) = dgi; %#ok<AGROW>
end

%% ------------------------------------------------------------------ %%
%  Structural check: G must be block diagonal (decentralised actuation)
%% ------------------------------------------------------------------ %%
% Cheap assertion, active only when the caller sets a debug flag, because it
% runs inside the integration loop.
if isfield(spec,'debug') && isfield(spec.debug,'checkBlockDiag') && spec.debug.checkBlockDiag
    mask = true(size(G));
    for i = 1:nA
        mask((i-1)*nx + (1:nx), (i-1)*nu + (1:nu)) = false;
    end
    assert(all(G(mask) == 0), ...
        'cpmg_dynamics_multiarea: G is not block diagonal; actuation is not decentralised');
end

end
