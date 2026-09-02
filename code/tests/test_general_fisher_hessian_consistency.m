function stats = test_general_fisher_hessian_consistency()
%TEST_GENERAL_FISHER_HESSIAN_CONSISTENCY Inline s, s', s'' Hessian test.

[problem, rho, direction] = makeProblem();
hessianDirection = src.solvers.FullHessianAction(rho, direction, problem);
tValues = 10 .^ (-(2:8));
errors = zeros(size(tValues));
for j = 1:numel(tValues)
    t = tValues(j);
    finiteDifference = (src.discretization.ps.Gradient( ...
        rho + t * direction, problem) - src.discretization.ps.Gradient( ...
        rho - t * direction, problem)) / (2 * t);
    difference = sqrt(problem.grid.h ...
        * sum((finiteDifference - hessianDirection) .^ 2));
    scale = max([1e-14, sqrt(problem.grid.h * sum(finiteDifference .^ 2)), ...
        sqrt(problem.grid.h * sum(hessianDirection .^ 2))]);
    errors(j) = difference / scale;
end
stats.best_relative_error = min(errors);
assert(stats.best_relative_error <= 1e-6, ...
    'General Fisher Hessian error %.3e.', stats.best_relative_error);
fprintf('test_general_fisher_hessian_consistency: %.3e\n', ...
    stats.best_relative_error);
end

function [problem, rho, direction] = makeProblem()
parameters = model.DefaultParameters1D();
parameters.L = 6;
parameters.N = 96;
grid = model.SetupGrid1D(parameters);
epsilon = parameters.epsilon;
amplitude = 5e-3;
width = 2e-2;
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = model.BuildPotential(grid, parameters.potential_label);
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = epsilon;
problem.regularization.transition_width = epsilon;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(x) ...
    x + epsilon + amplitude * (1 - exp(-x / width));
problem.fisher_regularization.ds_epsilon = @(x) ...
    1 + amplitude / width * exp(-x / width);
problem.fisher_regularization.d2s_epsilon = @(x) ...
    -amplitude / width ^ 2 * exp(-x / width);
problem.fisher_regularization.label = 'inline concave s_epsilon';
problem.potential_regularization = src.potential.MakeLinear();
rho = 0.04 + 0.02 * cos(pi * grid.x / grid.L + 0.2);
rho = rho / src.constraints.Mass(rho, grid.h);
direction = sin(3*pi*grid.x/grid.L) + 0.3*cos(5*pi*grid.x/grid.L);
direction = direction * (0.1 * min(rho) / max(abs(direction)));
end
