function stats = test_gradient_component_sum()
%TEST_GRADIENT_COMPONENT_SUM Diagnostic terms reproduce production G.

rng(41);
epsilon = 1e-3;
sigma = 1e-12;
parameters = model.DefaultParameters1D();
parameters.L = 8;
parameters.N = 128;
grid = model.SetupGrid1D(parameters);

problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = 0.5 * grid.x .^ 2;
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(rho) rho + epsilon;
problem.fisher_regularization.ds_epsilon = @(rho) ones(size(rho));
problem.fisher_regularization.d2s_epsilon = @(rho) zeros(size(rho));
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(rho) ...
    rho ./ hypot(rho, sigma);
problem.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
problem.potential_regularization.prox_type = 'generic_convex';

rho = 0.1 + rand(grid.N, 1);
rho = rho * (problem.mass / src.constraints.Mass(rho, grid.h));
components = src.diagnostics.GradientComponents(rho, problem);
fullGradient = src.discretization.ps.Gradient(rho, problem);
summed = components.fisher + components.potential ...
    + components.beta + components.delta;
stats.relative_error = norm(fullGradient - summed) ...
    / max(1, norm(fullGradient));
stats.total_field_error = norm(fullGradient - components.total) ...
    / max(1, norm(fullGradient));
stats.maximum_error = max(stats.relative_error, stats.total_field_error);
assert(stats.maximum_error <= 1e-13, ...
    'Gradient-component sum error is %.3e.', stats.maximum_error);
fprintf('test_gradient_component_sum: relative error %.3e\n', ...
    stats.maximum_error);
end
