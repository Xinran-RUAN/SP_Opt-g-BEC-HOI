function stats = test_gradient_consistency()
%TEST_GRADIENT_CONSISTENCY Centered directional finite differences.

problem = makeProblem();
x = problem.grid.x;
meanDensity = problem.mass / problem.grid.domain_length;
rho = meanDensity * (1 + 0.98 * cos(pi * (x + problem.grid.L) ...
    / problem.grid.L + 0.173));
rho = rho * (problem.mass / src.constraints.Mass(rho, problem.grid.h));

eta = sin(3 * pi * x / problem.grid.L) ...
    + 0.37 * cos(5 * pi * x / problem.grid.L + 0.21);
eta = eta - mean(eta);
eta = eta * (0.1 * min(rho) / (1e-2 * max(abs(eta))));
tValues = 10 .^ (-(2:8));
names = src.regularization.SupportedNames();
bestErrors = zeros(numel(names), 1);

for j = 1:numel(names)
    problem.regularization.name = names{j};
    gradient = src.discretization.ps.Gradient(rho, problem);
    exact = problem.grid.h * sum(gradient .* eta);
    errors = zeros(size(tValues));
    for ell = 1:numel(tValues)
        t = tValues(ell);
        plusEnergy = src.discretization.ps.Energy(rho + t * eta, problem);
        minusEnergy = src.discretization.ps.Energy(rho - t * eta, problem);
        finiteDifference = (plusEnergy - minusEnergy) / (2 * t);
        errors(ell) = abs(finiteDifference - exact) ...
            / max([abs(finiteDifference), abs(exact), 1e-12]);
    end
    bestErrors(j) = min(errors);
    assert(bestErrors(j) <= 1e-7, ...
        'Gradient error for %s is %.3e.', names{j}, bestErrors(j));
    fprintf('  gradient %-14s best relative error %.3e\n', ...
        names{j}, bestErrors(j));
end

stats.names = names;
stats.best_relative_error = bestErrors;
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
end
