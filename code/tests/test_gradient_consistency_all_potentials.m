function stats = test_gradient_consistency_all_potentials()
%TEST_GRADIENT_CONSISTENCY_ALL_POTENTIALS Target energy/gradient checks.

problem = makeProblem();
x = problem.grid.x;
rho = 0.04 + 0.02 * cos(pi * x / problem.grid.L + 0.17) ...
    + 0.005 * sin(3 * pi * x / problem.grid.L);
rho = rho * (problem.mass / src.constraints.Mass(rho, problem.grid.h));
direction = sin(3 * pi * x / problem.grid.L) ...
    + 0.37 * cos(5 * pi * x / problem.grid.L + 0.21);
direction = direction - mean(direction);
direction = direction * (0.1 * min(rho) / max(abs(direction)));
tValues = 10 .^ (-(2:8));
specifications = {
    struct('name', 'linear'), 'linear'
    struct('name', 'sqrt_same_scale'), 'legacy-p1'
    struct('name', 'sqrt_squared_scale'), 'legacy-p2'
    struct('name', 'sqrt_power', 'power', 1), 'sqrt-power-p1'
    struct('name', 'sqrt_power', 'power', 2), 'sqrt-power-p2'
    struct('name', 'sqrt_power', 'power', 3), 'sqrt-power-p3'
    struct('name', 'sqrt_power', 'power', 4), 'sqrt-power-p4'
};
labels = specifications(:, 2);
bestErrors = zeros(size(specifications, 1), 1);
for j = 1:size(specifications, 1)
    problem.potential_regularization = specifications{j, 1};
    gradient = src.discretization.ps.Gradient(rho, problem);
    exact = problem.grid.h * sum(gradient .* direction);
    errors = zeros(size(tValues));
    for q = 1:numel(tValues)
        t = tValues(q);
        plusEnergy = src.discretization.ps.Energy( ...
            rho + t * direction, problem);
        minusEnergy = src.discretization.ps.Energy( ...
            rho - t * direction, problem);
        finiteDifference = (plusEnergy - minusEnergy) / (2 * t);
        errors(q) = abs(finiteDifference - exact) ...
            / max([abs(finiteDifference), abs(exact), 1e-12]);
    end
    bestErrors(j) = min(errors);
    assert(bestErrors(j) <= 1e-7, ...
        'Potential-gradient error for %s is %.3e.', ...
        labels{j}, bestErrors(j));
    fprintf('  potential gradient %-20s best relative error %.3e\n', ...
        labels{j}, bestErrors(j));
end
stats.names = labels;
stats.best_relative_error = bestErrors;
stats.maximum_error = max(bestErrors);
end

function problem = makeProblem()
parameters = model.DefaultParameters1D();
parameters.L = 8;
parameters.N = 128;
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
problem.potential_regularization.name = 'linear';
end
