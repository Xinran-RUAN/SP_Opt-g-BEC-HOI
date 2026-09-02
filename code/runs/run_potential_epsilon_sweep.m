%RUN_POTENTIAL_EPSILON_SWEEP Manual P2 stiffness diagnostic; not automatic.
clearvars; clc;

epsilon_list = [1e-2; 5e-3; 2e-3; 1e-3];
potential_name = 'sqrt_squared_scale';
beta = 10;
delta = 10;
mass = 1;
L = 8;
N = 256;
main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
max_iter = 200000;
projection_name = 'semismooth';
save_result = true;
overwrite_existing = true;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
baseConfig = experiments.DefaultConfig();
baseConfig.parameters.beta = beta;
baseConfig.parameters.delta = delta;
baseConfig.parameters.mass = mass;
baseConfig.parameters.L = L;
baseConfig.parameters.N = N;
baseConfig.potential_regularization.name = potential_name;
baseConfig.entropy.enabled = false;
baseConfig.entropy.eta = 0;
baseConfig.solver.name = 'fista_cd';
baseConfig.solver.splitting = 'potential_prox';
baseConfig.solver.projection_name = projection_name;
baseConfig.solver.projection_tol = 1e-14;
baseConfig.solver.pg_tol = main_pg_switch_tol;
baseConfig.solver.final_pg_tol = final_pg_tol;
baseConfig.solver.certification_tol = final_pg_tol;
baseConfig.solver.residual_check_interval = 10;
baseConfig.solver.potential_prox.mass_tol = 1e-14;
baseConfig.solver.potential_prox.inner_tol = 1e-14;
baseConfig.solver.max_iter = max_iter;
baseConfig.solver.switch.enabled = true;
baseConfig.solver.switch.pg_entry_tol = main_pg_switch_tol;
baseConfig.solver.switch.max_main_iter = 20000;
baseConfig.solver.polish_mode = 'if_needed';
baseConfig.solver.polish.pg_tol = final_pg_tol;
baseConfig.solver.polish.gmres_tol = 5e-13;
baseConfig.solver.display = false;

grid = model.SetupGrid1D(baseConfig.parameters);
rhoInitial = exp(-grid.x .^ 2) / sqrt(pi);
rhoInitial = rhoInitial * (mass / src.constraints.Mass(rhoInitial, grid.h));
solutions = cell(numel(epsilon_list), 1);
records = repmat(struct('epsilon', NaN, 'sigma', NaN, ...
    'target_energy', NaN, 'baseline_energy', NaN, 'pg_residual', NaN, ...
    'iterations', NaN, 'elapsed_time', NaN, 'accepted_L_max', NaN, ...
    'mean_backtracks', NaN, 'zero_fraction_1e12', NaN, ...
    'fft_tail_ratio', NaN), numel(epsilon_list), 1);
for j = 1:numel(epsilon_list)
    config = baseConfig;
    config.parameters.epsilon = epsilon_list(j);
    config.regularization.name = 'shift_smooth';
    config.regularization.epsilon = epsilon_list(j);
    config.regularization.transition_width = epsilon_list(j);
    solutions{j} = experiments.SolveGroundState(config, rhoInitial);
    activeOptions.zero_tol_list = [0, 1e-14, 1e-12, 1e-10, 1e-8];
    active = src.diagnostics.ActiveSetDiagnostics( ...
        solutions{j}.rho, grid, activeOptions);
    records(j).epsilon = epsilon_list(j);
    records(j).sigma = solutions{j}.potential_regularization.sigma;
    records(j).target_energy = solutions{j}.target_energy;
    records(j).baseline_energy = solutions{j}.baseline_energy;
    records(j).pg_residual = solutions{j}.diagnostics.pg_residual;
    records(j).iterations = solutions{j}.diagnostics.iterations;
    records(j).elapsed_time = solutions{j}.diagnostics.elapsed_time;
    records(j).accepted_L_max = solutions{j}.diagnostics.accepted_L_max;
    records(j).mean_backtracks = solutions{j}.diagnostics.mean_backtracks;
    records(j).zero_fraction_1e12 = active.zero_fraction_by_tol(3);
    records(j).fft_tail_ratio = ...
        solutions{j}.diagnostics.fourier_tail.tail_ratio_quarter;
end
fprintf('\nP2 epsilon stiffness diagnostic\n');
fprintf('epsilon   sigma      target_E       baseline_E     PG       iter    time    max_L      mean_bt   zero1e-12   FFT_tail\n');
for j = 1:numel(records)
    r = records(j);
    fprintf('%.1e  %.1e  %.10e  %.10e  %.2e  %7d  %7.2f  %.2e  %.2f  %.3f  %.2e\n', ...
        r.epsilon, r.sigma, r.target_energy, r.baseline_energy, ...
        r.pg_residual, r.iterations, r.elapsed_time, r.accepted_L_max, ...
        r.mean_backtracks, r.zero_fraction_1e12, r.fft_tail_ratio);
end
archive.parameters = baseConfig.parameters;
archive.parameters.epsilon_list = epsilon_list;
archive.regularization.name = 'shift_smooth';
archive.potential_regularization.name = potential_name;
archive.grid = grid;
archive.solver = baseConfig.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = [records.target_energy].';
archive.target_energy = archive.energy;
archive.baseline_energy = [records.baseline_energy].';
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'potential_epsilon', ...
        sprintf('potential_epsilon_%s_N%d.mat', potential_name, N), ...
        overwrite_existing);
end
