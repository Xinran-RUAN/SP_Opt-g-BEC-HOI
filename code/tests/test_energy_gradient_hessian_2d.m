function stats = test_energy_gradient_hessian_2d()
%TEST_ENERGY_GRADIENT_HESSIAN_2D Directional consistency of the 2D core.

[problem, rho, eta] = makeProblem();
gradient = src.discretization.ps.Gradient(rho, problem);
exactEnergyDerivative = problem.grid.h * sum(gradient .* eta);
tValues = 10 .^ (-(2:8));
gradientErrors = nan(size(tValues));
hessianErrors = nan(size(tValues));
hessianEta = src.solvers.FullHessianAction(rho, eta, problem);
for j = 1:numel(tValues)
    t = tValues(j);
    energyDifference = (src.discretization.ps.Energy( ...
        rho + t * eta, problem) - src.discretization.ps.Energy( ...
        rho - t * eta, problem)) / (2 * t);
    gradientErrors(j) = abs(energyDifference - exactEnergyDerivative) ...
        / max([1e-14, abs(energyDifference), abs(exactEnergyDerivative)]);
    gradientDifference = (src.discretization.ps.Gradient( ...
        rho + t * eta, problem) - src.discretization.ps.Gradient( ...
        rho - t * eta, problem)) / (2 * t);
    difference = sqrt(problem.grid.h * sum( ...
        (gradientDifference - hessianEta) .^ 2));
    scale = max([1e-14, sqrt(problem.grid.h * sum(gradientDifference .^ 2)), ...
        sqrt(problem.grid.h * sum(hessianEta .^ 2))]);
    hessianErrors(j) = difference / scale;
end
stats.gradient_error = min(gradientErrors);
stats.hessian_error = min(hessianErrors);

rng(19);
u = randn(problem.grid.N, 1);
v = randn(problem.grid.N, 1);
Hu = src.solvers.FullHessianAction(rho, u, problem);
Hv = src.solvers.FullHessianAction(rho, v, problem);
uv = problem.grid.h * sum(u .* Hv);
vu = problem.grid.h * sum(Hu .* v);
stats.hessian_symmetry_error = abs(uv - vu) / max([1, abs(uv), abs(vu)]);
stats.rayleigh_quotient = problem.grid.h * sum(v .* Hv) ...
    / (problem.grid.h * sum(v .^ 2));

[P, preconditionerInfo] = src.solvers.BuildFDHessianPreconditioner(rho, problem);
stats.preconditioner_symmetry_error = preconditionerInfo.symmetry_error;
stats.preconditioner_min_diagonal = preconditionerInfo.min_diagonal;
setup.type = 'ict';
setup.droptol = 1e-3;
setup.diagcomp = 0;
R = ichol(P, setup); %#ok<NASGU>

assert(stats.gradient_error <= 1e-7, ...
    '2D energy-gradient error %.3e.', stats.gradient_error);
assert(stats.hessian_error <= 1e-6, ...
    '2D gradient-Hessian error %.3e.', stats.hessian_error);
assert(stats.hessian_symmetry_error <= 1e-11, ...
    '2D Hessian symmetry error %.3e.', stats.hessian_symmetry_error);
assert(stats.rayleigh_quotient > 0, '2D Hessian is not positive definite.');
assert(stats.preconditioner_symmetry_error <= 1e-13, ...
    '2D FD preconditioner symmetry error %.3e.', ...
    stats.preconditioner_symmetry_error);
fprintf(['test_energy_gradient_hessian_2d: gradient %.3e, Hessian %.3e, ' ...
    'symmetry %.3e, Rayleigh %.3e\n'], stats.gradient_error, ...
    stats.hessian_error, stats.hessian_symmetry_error, ...
    stats.rayleigh_quotient);
end

function [problem, rho, eta] = makeProblem()
parameters.L = 4;
parameters.Nx = 24;
parameters.Ny = 20;
grid = model.SetupGrid2D(parameters);
epsilon = 1e-2;
sigma = 1e-4;
problem.grid = grid;
problem.plan = src.discretization.ps.Plan2D(grid);
problem.V = 1 + 0.2 * cos(pi * grid.X / grid.L) ...
    + 0.1 * cos(2 * pi * grid.Y / grid.L);
problem.V = problem.V(:);
problem.beta = 10;
problem.delta = 10;
problem.mass = 1;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = epsilon;
problem.regularization.transition_width = epsilon;
problem.fisher_regularization.r_epsilon = @(x) x + epsilon;
problem.fisher_regularization.dr_epsilon = @(x) ones(size(x));
problem.fisher_regularization.d2r_epsilon = @(x) zeros(size(x));
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.label = 'r_epsilon(rho)=rho+epsilon';
problem.potential_regularization.p_sigma = @(x) ...
    x .^ 2 ./ (hypot(x, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(x) x ./ hypot(x, sigma);
problem.potential_regularization.d2p_sigma = @(x) ...
    (sigma ./ hypot(x, sigma)) .^ 2 ./ hypot(x, sigma);
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.label = 'sqrt potential smoothing';
problem.potential_regularization.name = 'custom_inline';
problem.potential_regularization.power = NaN;
problem.potential_regularization.prox_type = 'generic_convex';

rho = 1 + 0.12 * cos(pi * grid.X / grid.L) ...
    + 0.08 * sin(2 * pi * grid.Y / grid.L) ...
    + 0.03 * cos(pi * grid.X / grid.L) .* sin(pi * grid.Y / grid.L);
rho = rho(:);
rho = rho * (problem.mass / src.constraints.Mass(rho, grid.h));
eta = cos(pi * grid.X / grid.L) ...
    + 0.2 * cos(2 * pi * grid.Y / grid.L) ...
    + 0.1 * cos(pi * grid.X / grid.L) .* cos(pi * grid.Y / grid.L);
eta = eta(:) - mean(eta(:));
eta = eta * (0.1 * min(rho) / max(abs(eta)));
end
