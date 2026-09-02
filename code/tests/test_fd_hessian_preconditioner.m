function stats = test_fd_hessian_preconditioner()
%TEST_FD_HESSIAN_PRECONDITIONER SPD and iteration-count diagnostics.

problem = makeProblem(128);
x = problem.grid.x;
rho = exp(-x .^ 2) / sqrt(pi) + 1e-7;
rho = rho * (problem.mass / src.constraints.Mass(rho, problem.grid.h));
[P, info] = src.solvers.BuildFDHessianPreconditioner(rho, problem);
[R, cholFlag] = chol(P, 'lower');
assert(cholFlag == 0, 'Variable-FD preconditioner Cholesky failed.');
assert(info.symmetry_error <= 1e-14, ...
    'FD preconditioner symmetry error is %.3e.', info.symmetry_error);
rng(22);
minimumQuadraticForm = inf;
for j = 1:10
    v = randn(problem.grid.N, 1);
    minimumQuadraticForm = min(minimumQuadraticForm, v' * P * v);
end
assert(minimumQuadraticForm > 0, ...
    'FD preconditioner is not positive definite in sampled directions.');

rhs = randn(problem.grid.N, 1);
Hfun = @(v) src.solvers.FullHessianAction(rho, v, problem);
[~, noneFlag, noneRelres, noneIterations] = pcg( ...
    Hfun, rhs, 1e-8, 500);
Mfun = @(v) R' \ (R \ v);
[~, fdFlag, fdRelres, fdIterations] = pcg( ...
    Hfun, rhs, 1e-8, 500, Mfun);
assert(fdFlag == 0, ...
    'Variable-FD PCG failed with relres %.3e.', fdRelres);
if noneFlag == 0 && fdIterations >= 0.8 * noneIterations
    warning('test:fdPreconditioner:WeakImprovement', ...
        ['Variable-FD PCG iteration reduction was weak: ' ...
        'none=%d, FD=%d.'], noneIterations, fdIterations);
end

stats.cholesky_succeeded = cholFlag == 0;
stats.symmetry_error = info.symmetry_error;
stats.min_diagonal = info.min_diagonal;
stats.minimum_quadratic_form = minimumQuadraticForm;
stats.none_flag = noneFlag;
stats.none_relres = noneRelres;
stats.none_iterations = noneIterations;
stats.fd_flag = fdFlag;
stats.fd_relres = fdRelres;
stats.fd_iterations = fdIterations;
fprintf(['test_fd_hessian_preconditioner: chol=%d, symmetry %.3e, ' ...
    'PCG none/FD %d/%d\n'], cholFlag == 0, info.symmetry_error, ...
    noneIterations, fdIterations);
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
