function [P, info] = BuildFDHessianPreconditioner(rho, problem)
%BUILDFDHESSIANPRECONDITIONER Sparse variable-coefficient SPD surrogate.
%
% This matrix is used only as a PCG preconditioner. The target Hessian and
% the optimized Fourier pseudospectral discretization remain unchanged.

rho = rho(:);
N = problem.grid.N;
dimension = spatialDimension(problem.plan);
q = src.discretization.ps.SpatialGradient(rho, problem.plan);
[r, dr, d2r] = src.regularization.EvaluateFisher( ...
    rho, src.regularization.ResolveFisher(problem));
W = spdiags(1 ./ (4 * r), 0, N, N);
fisherRemainder = -sum(q .^ 2, 2) .* d2r ./ (8 * r .^ 2);

potential = problem.potential_regularization;
fisher = src.regularization.ResolveFisher(problem);
[~, ~, potentialSecondDerivative, ~] = src.potential.Evaluate( ...
    rho, potential, fisher.epsilon);
potentialCurvature = problem.V(:) .* potentialSecondDerivative;
diagonal = fisherRemainder + problem.beta + potentialCurvature;
if dimension == 1
    derivatives = {periodicForwardDifference(N, problem.grid.h)};
else
    Dx1 = periodicForwardDifference(problem.grid.Nx, problem.grid.hx);
    Dy1 = periodicForwardDifference(problem.grid.Ny, problem.grid.hy);
    derivatives = {kron(Dx1, speye(problem.grid.Ny)), ...
        kron(speye(problem.grid.Nx), Dy1)};
end

P = spdiags(diagonal, 0, N, N);
coefficient = dr ./ r;
for direction = 1:dimension
    Dfd = derivatives{direction};
    Bfd = Dfd - spdiags(q(:, direction) .* coefficient, 0, N, N);
    P = P + Bfd' * W * Bfd + problem.delta * (Dfd' * Dfd);
end
P = 0.5 * (P + P');

info.symmetry_error = norm(P - P', 'fro') / max(1, norm(P, 'fro'));
info.min_diagonal = min(diag(P));
info.preconditioner_shift = 0;
info.N = N;
info.dimension = dimension;
if dimension == 1
    info.factorization = 'chol';
else
    % A sparse incomplete Cholesky application keeps the two-dimensional
    % Newton preconditioner matrix-free with respect to the spectral target
    % Hessian and avoids a prohibitively dense/full factorization.
    info.factorization = 'ichol';
end
end

function D = periodicForwardDifference(N, h)
rows = [(1:N), (1:N)]';
columns = [(1:N), (2:N), 1]';
values = [-ones(N, 1); ones(N, 1)] / h;
D = sparse(rows, columns, values, N, N);
end

function dimension = spatialDimension(plan)
if isfield(plan, 'dimension')
    dimension = plan.dimension;
else
    dimension = 1;
end
end
