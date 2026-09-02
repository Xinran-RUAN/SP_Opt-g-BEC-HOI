function stats = test_inline_regularization_interface()
%TEST_INLINE_REGULARIZATION_INTERFACE Inline handles reproduce legacy P2.

[legacyProblem, rho, direction] = makeProblem();
epsilon = legacyProblem.regularization.epsilon;
sigma = 1e-6;
inlineProblem = legacyProblem;
inlineProblem.fisher_regularization.epsilon = epsilon;
inlineProblem.fisher_regularization.s_epsilon = @(x) x + epsilon;
inlineProblem.fisher_regularization.ds_epsilon = @(x) ones(size(x));
inlineProblem.fisher_regularization.d2s_epsilon = @(x) zeros(size(x));
inlineProblem.fisher_regularization.label = ...
    's_epsilon(rho) = rho + epsilon';
inlineProblem.potential_regularization.sigma = sigma;
inlineProblem.potential_regularization.p_sigma = @(x) ...
    x .^ 2 ./ (hypot(x, sigma) + sigma);
inlineProblem.potential_regularization.dp_sigma = @(x) ...
    x ./ hypot(x, sigma);
inlineProblem.potential_regularization.d2p_sigma = @(x) ...
    (sigma ./ hypot(x, sigma)) .^ 2 ./ hypot(x, sigma);
inlineProblem.potential_regularization.label = ...
    'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
inlineProblem.potential_regularization.prox_type = 'generic_convex';
inlineProblem.potential_regularization.name = 'inline_test';

legacyEnergy = src.discretization.ps.Energy(rho, legacyProblem);
inlineEnergy = src.discretization.ps.Energy(rho, inlineProblem);
legacyGradient = src.discretization.ps.Gradient(rho, legacyProblem);
inlineGradient = src.discretization.ps.Gradient(rho, inlineProblem);
[legacyProblem.potential_regularization, ~] = src.potential.Validate( ...
    legacyProblem.potential_regularization, epsilon, rho);
legacyHessian = src.solvers.FullHessianAction( ...
    rho, direction, legacyProblem);
inlineHessian = src.solvers.FullHessianAction( ...
    rho, direction, inlineProblem);

stats.energy_error = abs(legacyEnergy - inlineEnergy);
stats.gradient_error = norm(legacyGradient - inlineGradient, inf);
stats.hessian_error = norm(legacyHessian - inlineHessian, inf);
assert(stats.energy_error <= 1e-14);
assert(stats.gradient_error <= 1e-13);
assert(stats.hessian_error <= 1e-11);
fprintf('test_inline_regularization_interface: E %.3e, G %.3e, H %.3e\n', ...
    stats.energy_error, stats.gradient_error, stats.hessian_error);
end

function [problem, rho, direction] = makeProblem()
parameters = model.DefaultParameters1D();
parameters.L = 6;
parameters.N = 64;
grid = model.SetupGrid1D(parameters);
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = model.BuildPotential(grid, parameters.potential_label);
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = parameters.epsilon;
problem.regularization.transition_width = parameters.epsilon;
problem.potential_regularization.name = 'sqrt_squared_scale';
rho = 0.03 + exp(-grid.x .^ 2);
rho = rho / src.constraints.Mass(rho, grid.h);
direction = sin(pi * grid.x / grid.L) + 0.2 * cos(3*pi*grid.x/grid.L);
end
