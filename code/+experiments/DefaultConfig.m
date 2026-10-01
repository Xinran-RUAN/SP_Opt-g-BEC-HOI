function config = DefaultConfig()
%DEFAULTCONFIG Default configuration for one reliable 1D HOI solve.

config.parameters = model.DefaultParameters1D();
config.regularization.name = 'shift_smooth';
config.regularization.epsilon = config.parameters.epsilon;
config.regularization.transition_width = config.parameters.epsilon;

% Physical trapping potential V(x). Run files may replace this expression
% directly without modifying model.BuildPotential.
config.trapping_potential.V = @(x) 0.5 * x .^ 2;
config.trapping_potential.label = 'V(x) = x^2/2';

config.potential_regularization.name = 'linear';
config.potential_regularization.power = 0;
config.potential_regularization.sigma = 0;
config.potential_regularization.convexity_tol = 1e-14;

config.entropy.enabled = false;
config.entropy.eta = 0;
config.entropy.prox.mass_tol = 1e-14;
config.entropy.prox.lambda_max_iter = 50;
config.entropy.prox.inner_newton_max_iter = 30;
config.entropy.prox.inner_tol = 1e-14;

config.solver.name = 'spg';
% Backward-compatible global default. Formal potential studies override
% this with potential_prox while legacy_full_gradient remains available
% for solver A/B verification.
config.solver.splitting = 'legacy_full_gradient';
config.solver.projection_name = 'simplex';
config.solver.projection_tol = 1e-13;
config.solver.residual_step = 1;
config.solver.residual_check_interval = 1;
config.solver.active_tol = 1e-12;
config.solver.L0 = 1;
config.solver.backtrack_factor = 2;
config.solver.ista_L_decrease = 0.5;
config.solver.max_backtracks = 100;
config.solver.pg_tol = 1e-8;
config.solver.certification_tol = 1e-9;
config.solver.final_pg_tol = 1e-12;
config.solver.fallback_fista_iter = 0;
config.solver.energy_tol = 1e-14;
config.solver.max_iter = 200000;
% Optional canonical solver wall-clock limit.  Infinite by default, so
% existing production runs retain their iteration/stationarity behavior.
config.solver.time_limit = inf;
config.solver.a = 4;
config.solver.feasibility_tol = 1e-13;
config.solver.display = true;
config.solver.display_interval = 100;
config.solver.display_every = 200;
config.solver.switch.enabled = false;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = 1e-5;
config.solver.switch.max_main_iter = 20000;
config.solver.switch.forced_pg_tol = 1e-4;
% Production default: stop FISTA at the energy-window handoff and enter
% the configured KKT polish. Performance diagnostics may set this false
% to continue the identical FISTA iteration beyond the handoff.
config.solver.switch.stop_at_energy_handoff = true;
config.solver.history.capture_handoff_state = false;
config.solver.potential_prox.mass_tol = 1e-13;
config.solver.potential_prox.lambda_max_iter = 50;
config.solver.potential_prox.inner_tol = 1e-13;
config.solver.potential_prox.inner_max_iter = 30;
config.solver.potential_prox.bracket_max_iter = 60;
config.solver.spg_bb_type = 'bb1';
config.solver.spg_alpha_min = 1e-12;
config.solver.spg_alpha_max = 1e12;
config.solver.spg_alpha_reset = 1;
config.solver.spg_curvature_tol = 1e-18;
config.solver.spg_nonmonotone_M = 10;
config.solver.spg_c1 = 1e-4;
config.solver.spg_backtrack = 0.5;
config.solver.spg_min_lambda = 1e-14;
config.solver.spg_descent_tol = 1e-13;
config.solver.polish_mode = 'if_needed';
config.solver.polish.entry_pg_tol = 1e-5;
config.solver.polish.pg_tol = 1e-12;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.polish.pdas_max_iter = 50;
config.solver.polish.tau_active = 1;
config.solver.polish.active_step = 1;
config.solver.polish.active_tol = 1e-12;
config.solver.polish.gmres_tol = 1e-10;
config.solver.polish.gmres_maxit = 200;
config.solver.polish.gmres_failure_tol = 0.5;
config.solver.polish.pcg_tol_max = 1e-2;
config.solver.polish.pcg_tol_min = 1e-10;
config.solver.polish.pcg_forcing_factor = 0.1;
config.solver.polish.pcg_maxit = [];
config.solver.polish.fraction_to_boundary = 0.995;
config.solver.polish.residual_armijo = 1e-4;
config.solver.polish.stagnation_window = 5;
config.solver.polish.stagnation_rel_improvement = 1e-2;
config.solver.polish.acceptable_floor = 1e-10;
config.solver.polish.residual_c1 = 1e-4;
config.solver.polish.backtrack = 0.5;
config.solver.polish.min_step = 1e-12;
config.solver.polish.max_backtracks = 30;
config.solver.polish.mass_tol = 1e-12;
config.solver.polish.roundoff_tol = 1e-13;

config.output.show_plot = true;
config.output.save_result = true;
config.output.overwrite_existing = false;
end
