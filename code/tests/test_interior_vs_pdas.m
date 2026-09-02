function stats = test_interior_vs_pdas()
%TEST_INTERIOR_VS_PDAS Same-handoff minimizer comparison at N=64/128.

Nvalues = [64, 128];
records = repmat(struct(), numel(Nvalues), 1);
for j = 1:numel(Nvalues)
    [handoff, problem, solver] = handoffState(Nvalues(j));
    options = solver.polish;
    options.entry_pg_tol = 1e-4;
    options.pg_tol = 1e-10;
    options.display = false;
    options.max_iter = 20;
    options.pdas_max_iter = 50;

    interiorOptions = options;
    interiorOptions.linear_solver = 'interior_pcg_schur';
    interiorOptions.preconditioner = 'fd_variable';
    interior = src.solvers.PolishKKT( ...
        handoff, problem, interiorOptions);
    pdasOptions = options;
    pdasOptions.linear_solver = 'pdas_gmres';
    pdas = src.solvers.PolishKKT(handoff, problem, pdasOptions);

    stateDifference = sqrt(problem.grid.h ...
        * sum((interior.rho - pdas.rho) .^ 2));
    energyDifference = abs(interior.energy - pdas.energy);
    assert(~interior.switched_to_pdas, ...
        'Interior method unexpectedly used the PDAS fallback at N=%d.', ...
        Nvalues(j));
    assert(interior.pg_residual <= 1e-9, ...
        'Interior full PG is %.3e at N=%d.', ...
        interior.pg_residual, Nvalues(j));
    assert(pdas.pg_residual <= 1e-9, ...
        'PDAS full PG is %.3e at N=%d.', pdas.pg_residual, Nvalues(j));
    assert(stateDifference <= 1e-7, ...
        'Interior/PDAS state difference is %.3e at N=%d.', ...
        stateDifference, Nvalues(j));
    assert(energyDifference <= 1e-11, ...
        'Interior/PDAS energy difference is %.3e at N=%d.', ...
        energyDifference, Nvalues(j));

    records(j).N = Nvalues(j);
    records(j).state_L2_difference = stateDifference;
    records(j).energy_difference = energyDifference;
    records(j).interior_pg = interior.pg_residual;
    records(j).pdas_pg = pdas.pg_residual;
    records(j).interior_outer_iterations = interior.iterations;
    records(j).pdas_outer_iterations = pdas.iterations;
    records(j).pcg_z_iterations = interior.total_pcg_z_iterations;
    records(j).pcg_w_iterations = interior.total_pcg_w_iterations;
    records(j).gmres_iterations = pdas.total_gmres_iterations;
    fprintf(['  interior/PDAS N=%d: state %.3e, energy %.3e, ' ...
        'PG %.3e/%.3e\n'], Nvalues(j), stateDifference, ...
        energyDifference, interior.pg_residual, pdas.pg_residual);
end
stats.records = records;
end

function [rho, problem, solver] = handoffState(N)
config = experiments.DefaultConfig();
config.parameters.L = 8;
config.parameters.N = N;
config.parameters.epsilon = 1e-3;
config.parameters.beta = 10;
config.parameters.delta = 10;
config.regularization.name = 'shift_smooth';
config.regularization.epsilon = 1e-3;
config.regularization.transition_width = 1e-3;
config.potential_regularization.name = 'sqrt_squared_scale';
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.projection_tol = 1e-14;
config.solver.residual_check_interval = 10;
config.solver.max_iter = 20000;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = 1e-5;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'none';
config.solver.display = false;
grid = model.SetupGrid1D(config.parameters);
rho0 = exp(-grid.x .^ 2) / sqrt(pi);
rho0 = rho0 * (config.parameters.mass ...
    / src.constraints.Mass(rho0, grid.h));
handoff = experiments.SolveGroundState(config, rho0);
rho = handoff.rho;
problem = handoff.problem;
solver = handoff.solver;
end
