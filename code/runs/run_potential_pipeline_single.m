%RUN_POTENTIAL_PIPELINE_SINGLE Low-cost P2 FISTA-to-KKT validation.
clearvars; clc;

epsilon = 1e-3;
L = 8;
N = 128;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
config = experiments.DefaultConfig();
config.parameters.beta = 10;
config.parameters.delta = 10;
config.parameters.mass = 1;
config.parameters.L = L;
config.parameters.N = N;
config.parameters.epsilon = epsilon;
config.regularization.name = 'shift_smooth';
config.regularization.epsilon = epsilon;
config.regularization.transition_width = epsilon;
config.potential_regularization.name = 'sqrt_squared_scale';
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.projection_tol = 1e-14;
config.solver.residual_check_interval = 10;
config.solver.potential_prox.mass_tol = 1e-14;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.max_iter = 200000;
config.solver.final_pg_tol = 1e-12;
config.solver.certification_tol = 1e-12;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = 1e-5;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = 1e-12;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.display = true;
config.solver.display_every = 200;

grid = model.SetupGrid1D(config.parameters);
rhoInitial = exp(-grid.x .^ 2) / sqrt(pi);
rhoInitial = rhoInitial * (config.parameters.mass ...
    / src.constraints.Mass(rhoInitial, grid.h));
result = experiments.SolveGroundState(config, rhoInitial);
experiments.SaveResult(result, 'single', ...
    'potential_pipeline_P2_N128_eps_1em03.mat', true);
