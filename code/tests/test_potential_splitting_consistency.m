function stats = test_potential_splitting_consistency()
%TEST_POTENTIAL_SPLITTING_CONSISTENCY Verify f+potential equals target.

parameters = model.DefaultParameters1D();
parameters.L = 6;
parameters.N = 96;
parameters.epsilon = 1e-3;
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

x = grid.x;
rho = 0.05 + 0.01 * cos(pi * x / grid.L) ...
    + 0.004 * sin(4 * pi * x / grid.L);
rho = rho / src.constraints.Mass(rho, grid.h);
direction = sin(3 * pi * x / grid.L) ...
    + 0.2 * cos(7 * pi * x / grid.L);
direction = direction - mean(direction);
direction = direction * (0.1 * min(rho) / max(abs(direction)));
tValues = 10 .^ (-(2:8));
names = src.potential.SupportedNames();
maxDecompositionError = 0;
bestGradientErrors = zeros(numel(names), 1);
for j = 1:numel(names)
    [problem.potential_regularization, ~] = src.potential.Validate( ...
        struct('name', names{j}), parameters.epsilon);
    targetEnergy = src.discretization.ps.Energy(rho, problem);
    smoothEnergy = src.discretization.ps.SmoothEnergy(rho, problem);
    [potentialDensity, ~] = src.potential.Evaluate( ...
        rho, problem.potential_regularization, parameters.epsilon);
    potentialEnergy = grid.h * sum(problem.V .* potentialDensity);
    maxDecompositionError = max(maxDecompositionError, ...
        abs(targetEnergy - smoothEnergy - potentialEnergy));

    smoothGradient = src.discretization.ps.SmoothGradient(rho, problem);
    exact = grid.h * sum(smoothGradient .* direction);
    errors = zeros(size(tValues));
    for q = 1:numel(tValues)
        t = tValues(q);
        plus = src.discretization.ps.SmoothEnergy( ...
            rho + t * direction, problem);
        minus = src.discretization.ps.SmoothEnergy( ...
            rho - t * direction, problem);
        finiteDifference = (plus - minus) / (2 * t);
        errors(q) = abs(finiteDifference - exact) ...
            / max([1e-12, abs(finiteDifference), abs(exact)]);
    end
    bestGradientErrors(j) = min(errors);
end
assert(maxDecompositionError <= 5e-14, ...
    'Target splitting energy defect %.3e.', maxDecompositionError);
assert(max(bestGradientErrors) <= 1e-7, ...
    'Smooth gradient error %.3e.', max(bestGradientErrors));
stats.max_energy_decomposition_error = maxDecompositionError;
stats.best_gradient_error = bestGradientErrors;
fprintf(['test_potential_splitting_consistency: energy %.3e, ' ...
    'smooth gradient max %.3e\n'], ...
    maxDecompositionError, max(bestGradientErrors));
end
