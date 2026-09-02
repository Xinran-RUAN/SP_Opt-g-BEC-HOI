%RUN_ENTROPY_BIAS_RESOLUTION_TRADEOFF Fixed-N entropy overlap diagnostic.
clearvars; clc;

% -------------------------- editable settings --------------------------
regularization = 'shift_smooth';
epsilon = 1e-3;
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 256;
eta_list = [0; 2; 1; 5e-1; 2e-1; 1e-1; 5e-2; 2e-2; ...
    1e-2; 5e-3; 2e-3; 1e-3; 5e-4; 2e-4; 1e-4];
target_energy_bias = 1e-8;
min_transition_cells = 8;
layer_lower_relative_density = 1e-10;
layer_upper_relative_density = 1e-2;
pg_tol = 1e-8;
max_iter = 200000;
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
config.regularization.name = regularization;
config.regularization.epsilon = epsilon;
config.regularization.transition_width = epsilon;
config.solver.pg_tol = pg_tol;
config.solver.max_iter = max_iter;
config.solver.display = false;

config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'spg';
config.solver.polish_mode = 'if_needed';
baseline = experiments.SolveGroundState(config);
solutions = cell(numel(eta_list), 1);
solutions{1} = baseline;
warmStart = baseline.rho;
for j = 2:numel(eta_list)
    config.entropy.enabled = true;
    config.entropy.eta = eta_list(j);
    config.solver.name = 'fista_cd';
    config.solver.polish_mode = 'none';
    solutions{j} = experiments.SolveGroundState(config, warmStart);
    warmStart = solutions{j}.rho;
end

layerOptions.lower_relative_density = layer_lower_relative_density;
layerOptions.upper_relative_density = layer_upper_relative_density;
records = repmat(emptyRecord(), numel(eta_list), 1);
for j = 1:numel(eta_list)
    result = solutions{j};
    layer = src.diagnostics.EntropyLayerResolution( ...
        result.rho, result.grid, layerOptions);
    records(j).eta = eta_list(j);
    records(j).physical_energy = result.physical_energy;
    records(j).physical_energy_bias = ...
        result.physical_energy - baseline.physical_energy;
    records(j).density_L2_to_eta0 = sqrt(result.grid.h ...
        * sum((result.rho - baseline.rho) .^ 2));
    records(j).density_Linf_to_eta0 = ...
        max(abs(result.rho - baseline.rho));
    records(j).transition_cells = layer.transition_cells;
    records(j).transition_cells_left = layer.transition_cells_left;
    records(j).transition_cells_right = layer.transition_cells_right;
    records(j).fft_tail_ratio = ...
        result.diagnostics.fourier_tail.tail_ratio_quarter;
    records(j).prox_residual = result.diagnostics.pg_residual;
    records(j).min_density = min(result.rho);
    records(j).underflow_count = nnz(result.rho == 0);
    records(j).tail_mass = result.diagnostics.tail_mass;
    records(j).tail_max = result.diagnostics.tail_max;
    records(j).energy_bias_ok = j > 1 ...
        && records(j).physical_energy_bias >= -1e-10 ...
        && records(j).physical_energy_bias <= target_energy_bias;
    records(j).layer_resolved = j > 1 ...
        && records(j).transition_cells >= min_transition_cells;
    records(j).overlap = records(j).energy_bias_ok ...
        && records(j).layer_resolved ...
        && records(j).prox_residual <= 10 * pg_tol;
end

overlapIndices = find([records.overlap]);
overlap_exists = ~isempty(overlapIndices);
if overlap_exists
    overlap_eta = [records(overlapIndices).eta];
else
    overlap_eta = [];
end
if any([records(2:end).physical_energy_bias] < -1e-10)
    warning(['Significant negative physical energy bias detected; do not ' ...
        'interpret the overlap result until solver consistency is fixed.']);
    overlap_exists = false;
    overlap_eta = [];
end

fprintf('\nfixed-N entropy bias-resolution tradeoff: N=%d, h=%.3e\n', ...
    N, baseline.grid.h);
fprintf('eta       dEphys       L2-state     cells L/R   FFT-tail    prox-res     zeros   overlap\n');
for j = 1:numel(records)
    fprintf('%-8.1e  %+10.3e  %.3e   %3d/%-3d     %.3e   %.3e   %4d     %d\n', ...
        records(j).eta, records(j).physical_energy_bias, ...
        records(j).density_L2_to_eta0, ...
        records(j).transition_cells_left, ...
        records(j).transition_cells_right, records(j).fft_tail_ratio, ...
        records(j).prox_residual, records(j).underflow_count, ...
        records(j).overlap);
end
if overlap_exists
    fprintf('usable eta overlap exists:');
    fprintf(' %.3e', overlap_eta);
    fprintf('\n');
else
    fprintf(['no eta satisfies dEphys<=%.1e and transition_cells>=%d ' ...
        'with certified prox residual.\n'], ...
        target_energy_bias, min_transition_cells);
end

archive.parameters = config.parameters;
archive.parameters.eta_list = eta_list;
archive.parameters.target_energy_bias = target_energy_bias;
archive.parameters.min_transition_cells = min_transition_cells;
archive.parameters.layer_lower_relative_density = ...
    layer_lower_relative_density;
archive.parameters.layer_upper_relative_density = ...
    layer_upper_relative_density;
archive.regularization = baseline.regularization;
archive.grid = baseline.grid;
archive.solver = config.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = [records.physical_energy].';
archive.physical_energy = archive.energy;
archive.entropy_value = cellfun(@(s) s.entropy_value, solutions);
archive.augmented_energy = cellfun(@(s) s.augmented_energy, solutions);
archive.diagnostics.records = records;
archive.diagnostics.overlap_exists = overlap_exists;
archive.diagnostics.overlap_eta = overlap_eta;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'entropy_tradeoff', ...
        sprintf('entropy_tradeoff_N%d_eps_%s.mat', N, ...
        strrep(sprintf('%.0e', epsilon), '-', 'm')), overwrite_existing);
end
if show_plot
    makePlot(records, root, save_result, target_energy_bias, ...
        min_transition_cells);
end

function record = emptyRecord()
record.eta = NaN;
record.physical_energy = NaN;
record.physical_energy_bias = NaN;
record.density_L2_to_eta0 = NaN;
record.density_Linf_to_eta0 = NaN;
record.transition_cells = NaN;
record.transition_cells_left = NaN;
record.transition_cells_right = NaN;
record.fft_tail_ratio = NaN;
record.prox_residual = NaN;
record.min_density = NaN;
record.underflow_count = NaN;
record.tail_mass = NaN;
record.tail_max = NaN;
record.energy_bias_ok = false;
record.layer_resolved = false;
record.overlap = false;
end

function makePlot(records, root, saveFigure, targetBias, minimumCells)
positive = records(2:end);
eta = [positive.eta];
f = figure('Name', 'Entropy bias-resolution tradeoff');
tiledlayout(2, 2);
nexttile;
loglog(eta, abs([positive.physical_energy_bias]), 'o-'); hold on;
yline(targetBias, '--'); grid on; xlabel('eta'); ylabel('|dEphys|');
nexttile;
semilogx(eta, [positive.transition_cells], 's-'); hold on;
yline(minimumCells, '--'); grid on; xlabel('eta'); ylabel('transition cells');
nexttile;
loglog(eta, [positive.fft_tail_ratio], '^-');
grid on; xlabel('eta'); ylabel('FFT tail ratio');
nexttile;
loglog(eta, [positive.prox_residual], 'd-');
grid on; xlabel('eta'); ylabel('prox residual');
if saveFigure
    folder = fullfile(root, 'results', 'entropy_tradeoff');
    if ~isfolder(folder)
        mkdir(folder);
    end
    exportgraphics(f, fullfile(folder, 'entropy_bias_resolution_tradeoff.png'), ...
        'Resolution', 180);
end
end
