function stats = test_hessian_consistency()
%TEST_HESSIAN_CONSISTENCY Centered finite differences of the gradient.
%
% piecewise_c1 is interpreted with the selected generalized d2s used by
% the semismooth polish. The state/direction are kept away from its d2s
% jump, so this test does not claim global classical C2 regularity.

problem = makeProblem();
x = problem.grid.x;
meanDensity = problem.mass / problem.grid.domain_length;
phases = linspace(0.07, 1.3, 31);
bestGap = -inf;
rho = [];
for phase = phases
    candidate = meanDensity * (1 + 0.995 * cos( ...
        pi * (x + problem.grid.L) / problem.grid.L + phase));
    gap = min(abs(candidate - problem.regularization.transition_width));
    if gap > bestGap
        rho = candidate;
        bestGap = gap;
    end
end
rho = rho * (problem.mass / src.constraints.Mass(rho, problem.grid.h));

eta = sin(3 * pi * x / problem.grid.L + 0.13) ...
    + 0.29 * cos(5 * pi * x / problem.grid.L - 0.31);
distanceToBoundary = min(rho);
distanceToKink = min(abs(rho - problem.regularization.transition_width));
safeDistance = min(distanceToBoundary, distanceToKink);
eta = eta * (0.2 * safeDistance / (1e-2 * max(abs(eta))));
tValues = 10 .^ (-(2:8));
names = src.regularization.SupportedNames();
bestErrors = zeros(numel(names), 1);

for j = 1:numel(names)
    problem.regularization.name = names{j};
    hessianEta = src.solvers.HessianAction( ...
        rho, eta, problem, problem.plan, problem.regularization);
    errors = zeros(size(tValues));
    for ell = 1:numel(tValues)
        t = tValues(ell);
        plusGradient = src.discretization.ps.Gradient(rho + t * eta, problem);
        minusGradient = src.discretization.ps.Gradient(rho - t * eta, problem);
        finiteDifference = (plusGradient - minusGradient) / (2 * t);
        differenceNorm = discreteNorm(finiteDifference - hessianEta, problem.grid.h);
        scale = max([discreteNorm(finiteDifference, problem.grid.h), ...
            discreteNorm(hessianEta, problem.grid.h), 1e-14]);
        errors(ell) = differenceNorm / scale;
    end
    bestErrors(j) = min(errors);
    if strcmp(names{j}, 'shift_smooth')
        tolerance = 1e-7;
    elseif strcmp(names{j}, 'piecewise_c1')
        tolerance = 1e-5;
    else
        tolerance = 1e-6;
    end
    assert(bestErrors(j) <= tolerance, ...
        'Hessian error for %s is %.3e.', names{j}, bestErrors(j));
    fprintf('  hessian %-14s best relative error %.3e\n', ...
        names{j}, bestErrors(j));
end

stats.names = names;
stats.best_relative_error = bestErrors;
stats.distance_to_kink = distanceToKink;

problem.regularization.name = 'shift_smooth';
potentialSpecifications = {
    struct('name', 'linear'), 'linear'
    struct('name', 'sqrt_same_scale'), 'legacy-p1'
    struct('name', 'sqrt_squared_scale'), 'legacy-p2'
    struct('name', 'sqrt_power', 'power', 2), 'sqrt-power-p2'
    struct('name', 'sqrt_power', 'power', 3), 'sqrt-power-p3'
    struct('name', 'sqrt_power', 'power', 4), 'sqrt-power-p4'
};
potentialNames = potentialSpecifications(:, 2);
potentialErrors = zeros(size(potentialSpecifications, 1), 1);
for j = 1:size(potentialSpecifications, 1)
    [problem.potential_regularization, ~] = src.potential.Validate( ...
        potentialSpecifications{j, 1}, ...
        problem.regularization.epsilon);
    hessianEta = src.solvers.FullHessianAction(rho, eta, problem);
    errors = zeros(size(tValues));
    for ell = 1:numel(tValues)
        t = tValues(ell);
        plusGradient = src.discretization.ps.Gradient(rho + t * eta, problem);
        minusGradient = src.discretization.ps.Gradient(rho - t * eta, problem);
        finiteDifference = (plusGradient - minusGradient) / (2 * t);
        differenceNorm = discreteNorm( ...
            finiteDifference - hessianEta, problem.grid.h);
        scale = max([discreteNorm(finiteDifference, problem.grid.h), ...
            discreteNorm(hessianEta, problem.grid.h), 1e-14]);
        errors(ell) = differenceNorm / scale;
    end
    potentialErrors(j) = min(errors);
    assert(potentialErrors(j) <= 1e-6, ...
        'Full Hessian error for %s is %.3e.', ...
        potentialNames{j}, potentialErrors(j));
    fprintf('  full hessian %-18s best relative error %.3e\n', ...
        potentialNames{j}, potentialErrors(j));
end
stats.potential_names = potentialNames;
stats.potential_best_relative_error = potentialErrors;
end

function value = discreteNorm(vector, h)
value = sqrt(h * sum(vector(:) .^ 2));
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
