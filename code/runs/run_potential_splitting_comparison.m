%RUN_POTENTIAL_SPLITTING_COMPARISON P2 legacy/new FISTA splitting A/B.
clearvars; clc;

% -------------------------- fixed first study --------------------------
potential_name = 'sqrt_squared_scale';
epsilon = 1e-3;
beta = 10;
delta = 10;
mass = 1;
L = 8;
N_list = [64, 128];
splitting_names = {'legacy_full_gradient', 'potential_prox'};
pg_tol = 1e-10;
certification_tol = 1e-9;
max_iter = 200000;
projection_name = 'semismooth';
save_result = true;
overwrite_existing = true;
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
baseConfig = experiments.DefaultConfig();
baseConfig.parameters.beta = beta;
baseConfig.parameters.delta = delta;
baseConfig.parameters.mass = mass;
baseConfig.parameters.L = L;
baseConfig.parameters.epsilon = epsilon;
baseConfig.regularization.name = 'shift_smooth';
baseConfig.regularization.epsilon = epsilon;
baseConfig.regularization.transition_width = epsilon;
baseConfig.potential_regularization.name = potential_name;
baseConfig.entropy.enabled = false;
baseConfig.entropy.eta = 0;
baseConfig.solver.name = 'fista_cd';
baseConfig.solver.projection_name = projection_name;
baseConfig.solver.projection_tol = 1e-14;
baseConfig.solver.pg_tol = pg_tol;
baseConfig.solver.certification_tol = certification_tol;
baseConfig.solver.residual_check_interval = 10;
baseConfig.solver.max_iter = max_iter;
baseConfig.solver.polish_mode = 'none';
baseConfig.solver.display = false;
baseConfig.solver.potential_prox.mass_tol = 1e-14;
baseConfig.solver.potential_prox.inner_tol = 1e-14;

results = cell(numel(N_list), numel(splitting_names));
records = repmat(emptyRecord(), numel(N_list), numel(splitting_names));
stateDifference = nan(numel(N_list), 1);
energyDifference = nan(numel(N_list), 1);
bothCertified = false(numel(N_list), 1);
for n = 1:numel(N_list)
    config = baseConfig;
    config.parameters.N = N_list(n);
    grid = model.SetupGrid1D(config.parameters);
    rhoInitial = exp(-grid.x .^ 2) / sqrt(pi);
    rhoInitial = rhoInitial * (mass ...
        / src.constraints.Mass(rhoInitial, grid.h));
    fprintf('\nP2 splitting comparison: N=%d\n', N_list(n));
    for s = 1:numel(splitting_names)
        config.solver.splitting = splitting_names{s};
        fprintf('  solving %-22s ...\n', splitting_names{s});
        results{n, s} = experiments.SolveGroundState(config, rhoInitial);
        records(n, s) = makeRecord(results{n, s}, splitting_names{s});
        r = records(n, s);
        fprintf(['  done E=%.12e, iter=%d, time=%.2fs, ' ...
            'CompPG=%.3e, FullPG=%.3e, tau=%.3e\n'], ...
            r.target_energy, r.iterations, r.elapsed_time, ...
            r.composite_pg_residual, r.full_pg_residual, r.final_tau);
    end
    stateDifference(n) = sqrt(grid.h * sum( ...
        (results{n, 1}.rho - results{n, 2}.rho) .^ 2));
    energyDifference(n) = abs( ...
        results{n, 1}.target_energy - results{n, 2}.target_energy);
    bothCertified(n) = all([records(n, :).solver_ok]);
end

fprintf(['\nN    splitting               target_E          CompPG    FullPG    ' ...
    'iter      time      tau       Lmax      mean/max_bt  ' ...
    'mean/max_prox(lambda) active FFT_tail solverOK\n']);
for n = 1:numel(N_list)
    for s = 1:numel(splitting_names)
        r = records(n, s);
        fprintf(['%3d  %-22s  %.12e  %.2e  %.2e  %7d  %7.2f  ' ...
            '%.2e  %.2e  %.3f/%d       %.2f/%d             ' ...
            '%4d  %.2e  %d\n'], ...
            r.N, r.splitting, r.target_energy, ...
            r.composite_pg_residual, r.full_pg_residual, ...
            r.iterations, r.elapsed_time, r.final_tau, r.max_L, ...
            r.mean_backtracks, r.max_backtracks, ...
            r.mean_prox_lambda_iterations, ...
            r.max_prox_lambda_iterations, r.active_count, ...
            r.fft_tail_ratio, r.solver_ok);
    end
    fprintf(['     state_L2_difference=%.3e  energy_difference=%.3e  ' ...
        'both_certified=%d\n'], ...
        stateDifference(n), energyDifference(n), bothCertified(n));
end

archive.parameters = baseConfig.parameters;
archive.parameters.N_list = N_list;
archive.regularization.name = 'shift_smooth';
archive.potential_regularization.name = potential_name;
archive.grid = cellfun(@(r) r.grid, results, 'UniformOutput', false);
archive.solver = baseConfig.solver;
archive.rho = cellfun(@(r) r.rho, results, 'UniformOutput', false);
archive.energy = cellfun(@(r) r.target_energy, results);
archive.target_energy = archive.energy;
archive.baseline_energy = cellfun(@(r) r.baseline_energy, results);
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.diagnostics.state_L2_difference = stateDifference;
archive.diagnostics.energy_difference = energyDifference;
archive.diagnostics.both_certified = bothCertified;
archive.history = cellfun(@(r) r.history, results, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'potential_splitting', ...
        sprintf('potential_splitting_N64_128_eps_%s.mat', ...
        strrep(sprintf('%.0e', epsilon), '-', 'm')), overwrite_existing);
end

function record = makeRecord(result, splitting)
diagnostics = result.diagnostics;
record = emptyRecord();
record.N = result.grid.N;
record.splitting = splitting;
record.target_energy = result.target_energy;
record.composite_pg_residual = diagnostics.composite_pg_residual;
record.full_pg_residual = diagnostics.full_pg_residual;
record.iterations = diagnostics.iterations;
record.elapsed_time = diagnostics.elapsed_time;
record.final_tau = diagnostics.final_tau;
record.max_L = diagnostics.accepted_L_max;
record.mean_backtracks = diagnostics.mean_backtracks;
record.max_backtracks = diagnostics.max_backtracks;
record.mean_prox_lambda_iterations = ...
    diagnostics.mean_prox_lambda_iterations;
record.max_prox_lambda_iterations = ...
    diagnostics.max_prox_lambda_iterations;
record.mean_prox_inner_iterations = ...
    diagnostics.mean_prox_inner_iterations;
record.max_prox_inner_iterations = ...
    diagnostics.max_prox_inner_iterations;
record.active_count = nnz(result.rho == 0);
fourier = src.diagnostics.FourierTailDiagnostics(result.rho);
record.fft_tail_ratio = fourier.tail_ratio_quarter;
record.solver_ok = record.composite_pg_residual <= 10 * result.solver.pg_tol ...
    && record.full_pg_residual <= result.solver.certification_tol;
end

function record = emptyRecord()
record.N = NaN;
record.splitting = '';
record.target_energy = NaN;
record.composite_pg_residual = NaN;
record.full_pg_residual = NaN;
record.iterations = NaN;
record.elapsed_time = NaN;
record.final_tau = NaN;
record.max_L = NaN;
record.mean_backtracks = NaN;
record.max_backtracks = NaN;
record.mean_prox_lambda_iterations = NaN;
record.max_prox_lambda_iterations = NaN;
record.mean_prox_inner_iterations = NaN;
record.max_prox_inner_iterations = NaN;
record.active_count = NaN;
record.fft_tail_ratio = NaN;
record.solver_ok = false;
end
