function stats = test_entropy_gradient_hessian()
%TEST_ENTROPY_GRADIENT_HESSIAN Centered augmented-gradient consistency.

parameters = model.DefaultParameters1D();
parameters.L = 4;
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
eta = 0.07;
rho = 0.3 + 0.04 * cos(pi * grid.x / grid.L) ...
    + 0.02 * sin(3 * pi * grid.x / grid.L);
direction = sin(2 * pi * grid.x / grid.L + 0.2) ...
    + 0.3 * cos(5 * pi * grid.x / grid.L);
hessianDirection = src.solvers.HessianAction( ...
    rho, direction, problem, problem.plan, problem.regularization) ...
    + eta * direction ./ rho;
tValues = 10 .^ (-(2:7));
errors = zeros(size(tValues));
for j = 1:numel(tValues)
    t = tValues(j);
    plus = src.discretization.ps.Gradient(rho + t * direction, problem) ...
        + eta * src.entropy.Gradient(rho + t * direction);
    minus = src.discretization.ps.Gradient(rho - t * direction, problem) ...
        + eta * src.entropy.Gradient(rho - t * direction);
    finiteDifference = (plus - minus) / (2 * t);
    errors(j) = norm(finiteDifference - hessianDirection) ...
        / max(1, norm(hessianDirection));
end
stats.best_relative_error = min(errors);
assert(stats.best_relative_error <= 1e-7, ...
    'Entropy Hessian consistency error is %.3e.', stats.best_relative_error);
assert(max(abs(src.entropy.HessianDiagonal(rho) - 1 ./ rho)) == 0, ...
    'Entropy HessianDiagonal does not return 1/rho.');
fprintf('test_entropy_gradient_hessian: best relative %.3e\n', ...
    stats.best_relative_error);
end
