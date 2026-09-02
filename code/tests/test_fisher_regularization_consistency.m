function result = test_fisher_regularization_consistency()
%TEST_FISHER_REGULARIZATION_CONSISTENCY Energy/gradient/Hessian tests.

epsilon = 0.2;
parameters = model.DefaultParameters1D();
parameters.L = pi;
parameters.N = 64;
parameters.beta = 1.3;
parameters.delta = 0.7;
parameters.mass = 1;
grid = model.SetupGrid1D(parameters);
plan = src.discretization.ps.Plan1D(grid);
x = grid.x;

rho = epsilon * (1.15 + 0.65 * cos(x) + 0.20 * sin(2 * x));
nearMatch = abs(rho - epsilon) < 0.04 * epsilon;
rho(nearMatch) = rho(nearMatch) + 0.06 * epsilon;
assert(min(rho) > 0 && any(rho < epsilon) && any(rho > epsilon));
direction = sin(3 * x) + 0.3 * cos(5 * x);
direction = direction / max(abs(direction));

families = cell(3, 1);
families{1}.r_epsilon = @(z) z + epsilon;
families{1}.dr_epsilon = @(z) ones(size(z));
families{1}.d2r_epsilon = @(z) zeros(size(z));
families{1}.epsilon = epsilon;
families{1}.name = 'shift';
families{1}.label = 'r_epsilon(rho)=rho+epsilon';
families{2} = src.regularization.PiecewiseFisherCm(epsilon, 1);
families{3} = src.regularization.PiecewiseFisherCm(epsilon, 2);

problem.grid = grid;
problem.plan = plan;
problem.V = zeros(grid.N, 1);
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.potential_regularization = src.potential.MakeLinear();
problem.regularization.name = 'inline_handles';
steps = epsilon * 10 .^ (-(2:8));
gradientError = nan(3, 1);
hessianError = nan(3, 1);

for index = 1:3
    problem.fisher_regularization = families{index};
    gradient = src.discretization.ps.Gradient(rho, problem);
    hessianDirection = src.solvers.FullHessianAction( ...
        rho, direction, problem);
    directional = grid.h * sum(gradient .* direction);
    energyErrors = nan(size(steps));
    hessianErrors = nan(size(steps));
    for stepIndex = 1:numel(steps)
        step = steps(stepIndex);
        plus = rho + step * direction;
        minus = rho - step * direction;
        assert(min(plus) > 0 && min(minus) > 0);
        energyPlus = src.discretization.ps.Energy(plus, problem);
        energyMinus = src.discretization.ps.Energy(minus, problem);
        finiteDifference = (energyPlus - energyMinus) / (2 * step);
        energyErrors(stepIndex) = abs(finiteDifference - directional) ...
            / max([1, abs(finiteDifference), abs(directional)]);
        gradientPlus = src.discretization.ps.Gradient(plus, problem);
        gradientMinus = src.discretization.ps.Gradient(minus, problem);
        hessianDifference = (gradientPlus - gradientMinus) / (2 * step);
        hessianErrors(stepIndex) = norm( ...
            hessianDifference - hessianDirection) ...
            / max(1, norm(hessianDirection));
    end
    gradientError(index) = min(energyErrors);
    hessianError(index) = min(hessianErrors);
end

assert(max(gradientError) <= 1e-7, ...
    'Fisher energy/gradient consistency failed.');
assert(max(hessianError) <= 1e-6, ...
    'Fisher gradient/Hessian consistency failed.');
result.names = {'shift', 'piecewise_c1', 'piecewise_c2'};
result.gradient_relative_error = gradientError;
result.hessian_relative_error = hessianError;
fprintf(['test_fisher_regularization_consistency passed: ' ...
    'max grad %.3e, max Hessian %.3e.\n'], ...
    max(gradientError), max(hessianError));
end
