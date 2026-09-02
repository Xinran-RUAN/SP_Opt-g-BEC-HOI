function stats = test_interior_newton_pcg_schur()
%TEST_INTERIOR_NEWTON_PCG_SCHUR Compare PCG-Schur with a direct KKT solve.

problem = makeProblem(32);
x = problem.grid.x;
rho = problem.mass / problem.grid.domain_length ...
    * (1 + 0.15 * cos(pi * x / problem.grid.L) ...
    + 0.05 * sin(2 * pi * x / problem.grid.L));
rho = rho * (problem.mass / src.constraints.Mass(rho, problem.grid.h));
options = testOptions(problem.grid.N);
result = src.solvers.PolishInteriorNewtonPCG(rho, problem, options);
dPcg = result.first_newton_direction;
deltaLambdaPcg = result.first_delta_lambda;
assert(~isempty(dPcg), 'Interior PCG did not construct a Newton direction.');

gradient = src.discretization.ps.Gradient(rho, problem);
lambda = -(problem.grid.h * sum(gradient)) ...
    / (problem.grid.h * problem.grid.N);
stationarity = gradient + lambda;
massResidual = src.constraints.Mass(rho, problem.grid.h) - problem.mass;
N = problem.grid.N;
H = zeros(N, N);
for j = 1:N
    basis = zeros(N, 1);
    basis(j) = 1;
    H(:, j) = src.solvers.FullHessianAction(rho, basis, problem);
end
KKT = [H, ones(N, 1); problem.grid.h * ones(1, N), 0];
direct = KKT \ (-[stationarity; massResidual]);
dDirect = direct(1:N);
deltaLambdaDirect = direct(end);

directionError = norm(dPcg - dDirect) / max(norm(dDirect), realmin);
multiplierError = abs(deltaLambdaPcg - deltaLambdaDirect) ...
    / max(1, abs(deltaLambdaDirect));
massLinearizationError = abs(problem.grid.h * sum(dPcg) + massResidual);
assert(directionError <= 1e-8, ...
    'PCG-Schur direction relative error is %.3e.', directionError);
assert(multiplierError <= 1e-8, ...
    'PCG-Schur multiplier relative error is %.3e.', multiplierError);
assert(massLinearizationError <= 1e-10, ...
    'Schur mass linearization error is %.3e.', massLinearizationError);

rng(21);
maximumSymmetryError = 0;
minimumRayleigh = inf;
for trial = 1:8
    u = randn(N, 1);
    v = randn(N, 1);
    Hu = src.solvers.FullHessianAction(rho, u, problem);
    Hv = src.solvers.FullHessianAction(rho, v, problem);
    lhs = problem.grid.h * sum(u .* Hv);
    rhs = problem.grid.h * sum(Hu .* v);
    symmetryError = abs(lhs - rhs) / max([1, abs(lhs), abs(rhs)]);
    maximumSymmetryError = max(maximumSymmetryError, symmetryError);
    rayleigh = problem.grid.h * sum(v .* Hv) ...
        / (problem.grid.h * sum(v .^ 2));
    minimumRayleigh = min(minimumRayleigh, rayleigh);
end
assert(maximumSymmetryError <= 1e-11, ...
    'Hessian symmetry error is %.3e.', maximumSymmetryError);
assert(minimumRayleigh > 0, ...
    'Hessian minimum sampled Rayleigh quotient is %.3e.', minimumRayleigh);

stats.direction_relative_error = directionError;
stats.multiplier_relative_error = multiplierError;
stats.mass_linearization_residual = massLinearizationError;
stats.max_symmetry_error = maximumSymmetryError;
stats.min_rayleigh_quotient = minimumRayleigh;
fprintf(['test_interior_newton_pcg_schur: direction %.3e, ' ...
    'multiplier %.3e, mass %.3e, symmetry %.3e, Rayleigh %.3e\n'], ...
    directionError, multiplierError, massLinearizationError, ...
    maximumSymmetryError, minimumRayleigh);
end

function options = testOptions(N)
options.pg_tol = 1e-30;
options.max_iter = 1;
options.active_step = 1e-8;
options.active_tol = 0;
options.projection_tol = 1e-14;
options.residual_step = 1;
options.pcg_tol_max = 1e-12;
options.pcg_tol_min = 1e-12;
options.pcg_forcing_factor = 0.1;
options.pcg_maxit = max(100, N);
options.preconditioner = 'fd_variable';
options.display = false;
end

function problem = makeProblem(N)
parameters = model.DefaultParameters1D();
parameters.L = 8;
parameters.N = N;
parameters.epsilon = 1e-3;
grid = model.SetupGrid1D(parameters);
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = model.BuildPotential(grid, parameters.potential_label);
problem.beta = 10;
problem.delta = 10;
problem.mass = 1;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = 1e-3;
problem.regularization.transition_width = 1e-3;
[problem.potential_regularization, ~] = src.potential.Validate( ...
    struct('name', 'sqrt_squared_scale'), 1e-3);
end
