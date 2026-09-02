%RUN_POTENTIAL_REGULARIZATION_COMPARISON Fixed-grid P0/P1/P2 comparison.
clearvars; clc;

% -------------------------- fixed first study --------------------------
epsilon = 1e-3;
beta = 10;
delta = 10;
mass = 1;
L = 8;
N = 256;
potential_names = {
    'linear'
    'sqrt_same_scale'
    'sqrt_squared_scale'
};
solver_name = 'fista_cd';
projection_name = 'semismooth';
main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
max_iter = 200000;
target_baseline_energy_bias = 1e-8; % display reference only
show_plot = true;
save_result = true;
overwrite_existing = true;
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
config = experiments.DefaultConfig();
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.mass = mass;
config.parameters.L = L;
config.parameters.N = N;
config.parameters.epsilon = epsilon;
config.regularization.name = 'shift_smooth';
config.regularization.epsilon = epsilon;
config.regularization.transition_width = epsilon;
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = solver_name;
config.solver.splitting = 'potential_prox';
config.solver.projection_name = projection_name;
config.solver.projection_tol = 1e-14;
config.solver.pg_tol = main_pg_switch_tol;
config.solver.final_pg_tol = final_pg_tol;
config.solver.certification_tol = final_pg_tol;
config.solver.residual_check_interval = 10;
config.solver.potential_prox.mass_tol = 1e-14;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.max_iter = max_iter;
config.solver.switch.enabled = true;
config.solver.switch.pg_entry_tol = main_pg_switch_tol;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = final_pg_tol;
config.solver.polish.gmres_tol = 5e-13;
config.solver.display = false;

grid = model.SetupGrid1D(config.parameters);
rhoInitial = exp(-grid.x .^ 2) / sqrt(pi);
rhoInitial = rhoInitial * (mass / src.constraints.Mass(rhoInitial, grid.h));
if min(rhoInitial) <= 0
    error('The common Gaussian initial state must be strictly positive.');
end

solutions = cell(numel(potential_names), 1);
for j = 1:numel(potential_names)
    variantConfig = config;
    variantConfig.potential_regularization.name = potential_names{j};
    solutions{j} = experiments.SolveGroundState(variantConfig, rhoInitial);
end
baselineResult = solutions{1};
records = repmat(src.diagnostics.PotentialStudyDiagnostics( ...
    solutions{1}, baselineResult), numel(potential_names), 1);
for j = 2:numel(potential_names)
    records(j) = src.diagnostics.PotentialStudyDiagnostics( ...
        solutions{j}, baselineResult);
end
if any([records.baseline_energy_bias] < -1e-9)
    warning('Baseline or modified solve is not sufficiently converged.');
end

fprintf('\npotential regularization comparison: epsilon=%.1e, L=%g, N=%d\n', ...
    epsilon, L, N);
fprintf(['variant                 sigma       target_E          baseline_E        ' ...
    'base_bias    PG         iter     time    zeros  near1e-12  tail_mass   FFT_tail    max_L      mean/max bt  restart\n']);
for j = 1:numel(records)
    r = records(j);
    fprintf(['%-22s  %.1e  %.12e  %.12e  %+9.2e  %.2e  %7d  %7.2f  ' ...
        '%5d   %.3f      %.2e   %.2e   %.2e   %.2f/%d      %d\n'], ...
        r.variant, r.sigma, r.target_energy, r.baseline_energy, ...
        r.baseline_energy_bias, r.pg_residual, r.iterations, ...
        r.elapsed_time, r.exact_zero_count, r.zero_fraction_1e12, ...
        r.tail_mass, r.fft_tail_ratio_quarter, r.accepted_L_max, ...
        r.mean_backtracks, r.max_backtracks, r.restart_count);
end

p0 = records(1);
p2 = records(3);
activeImproved = p2.zero_fraction_1e12 < 0.5 * p0.zero_fraction_1e12;
fourierImproved = p2.fft_tail_ratio_quarter ...
    < 0.5 * p0.fft_tail_ratio_quarter;
farFieldSmall = p2.boundary_strip_max <= 1e-12 ...
    && p2.boundary_strip_max_abs_drho <= 1e-10;
stationaritySmall = p2.pg_residual <= 10 * final_pg_tol;
fprintf(['P2 evidence: active_improved=%d, Fourier_improved=%d, ' ...
    'far_field_small=%d, PG_certified=%d.\n'], ...
    activeImproved, fourierImproved, farFieldSmall, stationaritySmall);
fprintf('P2 baseline bias %.3e; display reference %.3e.\n', ...
    p2.baseline_energy_bias, target_baseline_energy_bias);
if activeImproved && fourierImproved && farFieldSmall && stationaritySmall
    candidate_status = ['structural diagnostics improved; baseline-bias ' ...
        'acceptance is still required before mesh refinement'];
else
    candidate_status = 'inconclusive';
end
fprintf('P2 candidate status: %s\n', candidate_status);

archive.parameters = config.parameters;
archive.parameters.target_baseline_energy_bias = ...
    target_baseline_energy_bias;
archive.regularization = baselineResult.regularization;
archive.potential_regularization.names = potential_names;
archive.grid = grid;
archive.solver = config.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = [records.target_energy].';
archive.target_energy = archive.energy;
archive.baseline_energy = [records.baseline_energy].';
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.diagnostics.candidate_status = candidate_status;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'potential_comparison', ...
        sprintf('potential_comparison_N%d_eps_%s.mat', N, ...
        strrep(sprintf('%.0e', epsilon), '-', 'm')), overwrite_existing);
end
if show_plot && usejava('desktop')
    makePlots(solutions, records, root, save_result, ...
        target_baseline_energy_bias);
elseif show_plot
    warning(['Plots were skipped because MATLAB has no interactive desktop; ' ...
        'the four-figure path remains enabled for an interactive run.']);
end

function makePlots(solutions, records, root, saveFigures, biasReference)
labels = {records.variant};
colors = lines(numel(records));
folder = fullfile(root, 'results', 'potential_comparison');
if saveFigures && ~isfolder(folder)
    mkdir(folder);
end

f1 = figure('Name', 'Potential regularization densities'); hold on;
for j = 1:numel(records)
    plot(solutions{j}.grid.x, solutions{j}.rho, ...
        'Color', colors(j, :), 'LineWidth', 1.2);
end
xlabel('x'); ylabel('rho'); grid on; legend(labels, 'Location', 'best');
saveFigure(f1, folder, '01_density.png', saveFigures);

f2 = figure('Name', 'Potential regularization semilog tails'); hold on;
for j = 1:numel(records)
    relativeDensity = solutions{j}.rho / max(solutions{j}.rho);
    rhoPlot = max(relativeDensity, 1e-300); % plotting only
    semilogy(solutions{j}.grid.x, rhoPlot, ...
        'Color', colors(j, :), 'LineWidth', 1.2);
end
xlabel('x'); ylabel('rho/max(rho)'); grid on;
legend(labels, 'Location', 'best');
saveFigure(f2, folder, '02_semilog_density.png', saveFigures);

f3 = figure('Name', 'Potential regularization Fourier decay'); hold on;
for j = 1:numel(records)
    semilogy(abs(records(j).fourier_modes), records(j).rho_hat_abs, '.', ...
        'Color', colors(j, :));
end
set(gca, 'YScale', 'log'); xlabel('|k|'); ylabel('abs(rho hat)');
grid on; legend(labels, 'Location', 'best');
saveFigure(f3, folder, '03_fourier_decay.png', saveFigures);

f4 = figure('Name', 'Potential bias and FISTA cost');
tiledlayout(1, 2);
nexttile;
semilogy(1:numel(records), abs([records.baseline_energy_bias]), 'o-');
yline(biasReference, '--'); xticks(1:numel(records)); xticklabels(labels);
xtickangle(20); ylabel('|baseline energy bias|'); grid on;
nexttile;
yyaxis left; bar(1:numel(records), [records.iterations]);
ylabel('iterations');
yyaxis right; plot(1:numel(records), [records.elapsed_time], 'o-');
ylabel('elapsed time'); xticks(1:numel(records)); xticklabels(labels);
xtickangle(20); grid on;
saveFigure(f4, folder, '04_bias_solver_cost.png', saveFigures);
end

function saveFigure(handle, folder, name, enabled)
if enabled
    exportgraphics(handle, fullfile(folder, name), 'Resolution', 180);
end
end
