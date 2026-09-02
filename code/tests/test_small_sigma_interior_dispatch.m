function stats = test_small_sigma_interior_dispatch()
%TEST_SMALL_SIGMA_INTERIOR_DISPATCH Projected active is diagnostic only.

sigmaList = [1e-9, 1e-12];
stats.projected_active_count = zeros(size(sigmaList));
for j = 1:numel(sigmaList)
    [problem, rho, options] = makeCase(sigmaList(j));
    diagnostic = src.solvers.InteriorKKTResidual( ...
        rho, problem, options);
    assert(all(rho > 0) && diagnostic.exact_zero_count == 0);
    assert(diagnostic.projected_active_count > 0, ...
        'Synthetic handoff state did not exercise projected-active logic.');
    result = src.solvers.PolishKKT(rho, problem, options);
    assert(strcmp(result.polish_solver, 'interior_newton_pcg'), ...
        'Small-sigma strictly positive state incorrectly selected PDAS.');
    assert(~result.switched_to_pdas);
    stats.projected_active_count(j) = diagnostic.projected_active_count;
end
stats.sigma = sigmaList;
fprintf('test_small_sigma_interior_dispatch: projAct %d/%d, interior PASS\n', ...
    stats.projected_active_count(1), stats.projected_active_count(2));
end

function [problem, rho, options] = makeCase(sigma)
parameters = model.DefaultParameters1D();
parameters.L = 8;
parameters.N = 64;
grid = model.SetupGrid1D(parameters);
epsilon = 1e-3;
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = model.BuildPotential(grid, parameters.potential_label);
problem.beta = 10;
problem.delta = 10;
problem.mass = 1;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = epsilon;
problem.regularization.transition_width = epsilon;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(x) x + epsilon;
problem.fisher_regularization.ds_epsilon = @(x) ones(size(x));
problem.fisher_regularization.d2s_epsilon = @(x) zeros(size(x));
problem.fisher_regularization.label = 's_epsilon(rho) = rho + epsilon';
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.p_sigma = @(x) ...
    x .^ 2 ./ (hypot(x, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(x) x ./ hypot(x, sigma);
problem.potential_regularization.d2p_sigma = @(x) ...
    (sigma ./ hypot(x, sigma)) .^ 2 ./ hypot(x, sigma);
problem.potential_regularization.label = 'inline small-sigma potential';
problem.potential_regularization.prox_type = 'generic_convex';
problem.potential_regularization.name = 'inline_test';
rho = exp(-grid.x .^ 2) + sigma;
rho = rho / src.constraints.Mass(rho, grid.h);
options.residual_step = 1;
options.projection_tol = 1e-14;
options.active_tol = 1e-12;
options.linear_solver = 'interior_pcg_schur';
options.preconditioner = 'fd_variable';
options.allow_pdas_fallback = true;
options.max_iter = 0; % dispatch-only unit test
options.display = false;
options.pg_tol = 1e-12;
end
