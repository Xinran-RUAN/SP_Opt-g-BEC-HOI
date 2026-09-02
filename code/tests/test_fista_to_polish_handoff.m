function stats = test_fista_to_polish_handoff()
%TEST_FISTA_TO_POLISH_HANDOFF Energy-window handoff on the stiff P2 case.

config = experiments.DefaultConfig();
config.parameters.L = 8;
config.parameters.N = 64;
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
config.solver.residual_step = 1;
config.solver.residual_check_interval = 10;
config.solver.max_iter = 20000;
config.solver.pg_tol = 1e-10;
config.solver.final_pg_tol = 1e-12;
config.solver.certification_tol = 1e-9;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = 1e-5;
config.solver.switch.max_main_iter = 20000;
config.solver.switch.forced_pg_tol = 1e-4;
config.solver.history.capture_handoff_state = true;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = 1e-12;
config.solver.polish.max_iter = 50;
config.solver.display = false;
config.solver.potential_prox.mass_tol = 1e-14;
config.solver.potential_prox.inner_tol = 1e-14;

grid = model.SetupGrid1D(config.parameters);
rhoInitial = exp(-grid.x .^ 2) / sqrt(pi);
rhoInitial = rhoInitial * (config.parameters.mass ...
    / src.constraints.Mass(rhoInitial, grid.h));

pipeline = experiments.SolveGroundState(config, rhoInitial);
d = pipeline.diagnostics;
assert(strcmp(d.main_stop_reason, 'switch_to_kkt_polish'), ...
    'FISTA did not stop at the energy-window handoff.');
assert(d.main_pg_residual <= config.solver.switch.pg_entry_tol, ...
    'Handoff full PG %.3e exceeds its entry tolerance.', ...
    d.main_pg_residual);
assert(d.polish_entered, 'KKT polish was not entered.');
fprintf(['  handoff: stop=%s, main=%d, mainPG=%.3e, ' ...
    'polish=%s/%d, finalPG=%.3e, finalKKT=%.3e, dE=%.3e\n'], ...
    d.main_stop_reason, d.main_iterations, d.main_pg_residual, ...
    d.polish_status, d.polish_iterations, d.final_pg_residual, ...
    d.final_kkt_residual, d.polish_energy_change);
if d.polish_failed
    fprintf('  polish failure: %s\n', d.polish_failure_message);
    disp(struct2table(pipeline.history.polish));
end
assert(d.final_pg_residual < 1e-2 * d.main_pg_residual, ...
    'KKT polish did not materially reduce full PG.');
assert(d.final_pg_residual <= 10 * config.solver.polish.pg_tol, ...
    'Final full PG %.3e misses the polish tolerance.', ...
    d.final_pg_residual);
assert(d.polish_energy_change <= 1e-10 * max(1, abs(d.main_energy)), ...
    'KKT polish produced a material target-energy increase.');
assert(d.mass_error <= 1e-12, 'Pipeline mass error is %.3e.', d.mass_error);
assert(min(pipeline.rho) >= -1e-13, ...
    'Pipeline positivity error is %.3e.', min(pipeline.rho));

pureConfig = config;
pureConfig.solver.polish_mode = 'none';
pureConfig.solver.switch.stop_at_energy_handoff = false;
pureConfig.solver.final_pg_tol = 1e-16;
pureConfig.solver.max_iter = d.main_iterations + 20;
pure = experiments.SolveGroundState(pureConfig, rhoInitial);
assert(pure.diagnostics.handoff_iteration == d.main_iterations, ...
    'Diagnostic FISTA changed the energy-window handoff iteration.');
assert(pure.diagnostics.main_iterations > d.main_iterations, ...
    'Diagnostic FISTA did not continue beyond the energy handoff.');
handoffDifference = pure.diagnostics.handoff_rho ...
    - pipeline.diagnostics.handoff_rho;
assert(sqrt(grid.h * sum(handoffDifference .^ 2)) <= 1e-13, ...
    'Diagnostic and production handoff states differ.');
assert(max(abs(pure.history.main.augmented_energy(1:d.main_iterations) ...
    - pipeline.history.main.augmented_energy)) <= 1e-13, ...
    'Diagnostic FISTA changed the pre-handoff energy history.');
assert(max(abs(pure.history.main.full_pg_residual(1:d.main_iterations) ...
    - pipeline.history.main.full_pg_residual)) <= 1e-13, ...
    'Diagnostic FISTA changed the pre-handoff PG history.');
assert(d.total_iterations < pure.diagnostics.iterations, ...
    'Pipeline iteration count was not lower than pure FISTA.');

stats.handoff_iteration = d.main_iterations;
stats.handoff_energy = d.main_energy;
stats.window_energy_span = d.main_energy_window_span;
stats.handoff_pg = d.main_pg_residual;
stats.polish_iterations = d.polish_iterations;
stats.final_pg = d.final_pg_residual;
stats.final_kkt = d.final_kkt_residual;
stats.polish_energy_change = d.polish_energy_change;
stats.polish_state_change = d.polish_state_change;
stats.pipeline_iterations = d.total_iterations;
stats.pipeline_time = d.total_elapsed_time;
stats.pure_fista_iterations = pure.diagnostics.iterations;
stats.pure_fista_time = pure.diagnostics.total_elapsed_time;
fprintf(['  FISTA->polish: main=%d, polish=%d, PG %.3e -> %.3e, ' ...
    'pure FISTA=%d iterations\n'], d.main_iterations, ...
    d.polish_iterations, d.main_pg_residual, d.final_pg_residual, ...
    pure.diagnostics.iterations);
end
