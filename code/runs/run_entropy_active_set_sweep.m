%RUN_ENTROPY_ACTIVE_SET_SWEEP Diagnose entropy, active sets, and FFT tails.
clearvars; clc;

% -------------------------- fixed first study --------------------------
regularization = 'shift_smooth';
epsilon = 1e-3;
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 256;
eta_list = [0; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6];

pg_tol = 1e-8;
max_iter = 200000;
polish_mode = 'if_needed';
polish_kkt_tol = 1e-10;
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
config.solver.polish_mode = polish_mode;
config.solver.polish.pg_tol = polish_kkt_tol;
config.solver.display = false;

% The disabled and enabled-with-eta-zero solves start identically. Their
% direct comparison certifies that eta=0 never enters the entropy prox.
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'spg';
baselineDisabled = experiments.SolveGroundState(config);
config.entropy.enabled = true;
baselineEtaZero = experiments.SolveGroundState(config);
eta_zero_bitwise_state = isequal(baselineDisabled.rho, baselineEtaZero.rho);
eta_zero_bitwise_energy = isequal( ...
    baselineDisabled.physical_energy, baselineEtaZero.physical_energy);

numberOfEta = numel(eta_list);
solutions = cell(numberOfEta, 1);
solutions{1} = baselineEtaZero;
warmStart = baselineEtaZero.rho;
for etaIndex = 2:numberOfEta
    config.entropy.enabled = true;
    config.entropy.eta = eta_list(etaIndex);
    config.solver.name = 'fista_cd';
    solutions{etaIndex} = experiments.SolveGroundState(config, warmStart);
    warmStart = solutions{etaIndex}.rho;
end

% A second eta=0 solve comes from the final continuation state.
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'spg';
etaZeroFromContinuation = experiments.SolveGroundState(config, warmStart);
eta_zero_continuation_L2 = sqrt(baselineEtaZero.grid.h ...
    * sum((etaZeroFromContinuation.rho - baselineEtaZero.rho) .^ 2));
eta_zero_continuation_Linf = max(abs( ...
    etaZeroFromContinuation.rho - baselineEtaZero.rho));
eta_zero_continuation_energy_difference = ...
    etaZeroFromContinuation.physical_energy - baselineEtaZero.physical_energy;

records = repmat(emptyRecord(), numberOfEta, 1);
rhoZero = baselineEtaZero.rho;
Ezero = baselineEtaZero.physical_energy;
for etaIndex = 1:numberOfEta
    result = solutions{etaIndex};
    active = result.diagnostics.active_set;
    fourier = result.diagnostics.fourier_tail;
    [~, nearIndex] = min(abs(active.zero_tol_list - 1e-12));
    records(etaIndex).eta = eta_list(etaIndex);
    records(etaIndex).physical_energy = result.physical_energy;
    records(etaIndex).entropy_value = result.entropy_value;
    records(etaIndex).augmented_energy = result.augmented_energy;
    records(etaIndex).physical_energy_bias = result.physical_energy - Ezero;
    records(etaIndex).density_L2_to_eta0 = sqrt(result.grid.h ...
        * sum((result.rho - rhoZero) .^ 2));
    records(etaIndex).density_Linf_to_eta0 = max(abs(result.rho - rhoZero));
    records(etaIndex).min_density = min(result.rho);
    records(etaIndex).min_positive_density = active.min_positive_density;
    records(etaIndex).exact_zero_count = active.exact_zero_count;
    records(etaIndex).zero_fraction = active.zero_fraction_by_tol(nearIndex);
    records(etaIndex).zero_tol_list = active.zero_tol_list;
    records(etaIndex).zero_count_by_tol = active.zero_count_by_tol;
    records(etaIndex).zero_fraction_by_tol = active.zero_fraction_by_tol;
    records(etaIndex).support_left_x = active.support_left_x;
    records(etaIndex).support_right_x = active.support_right_x;
    records(etaIndex).support_width = active.support_width;
    records(etaIndex).fft_tail_ratio = fourier.tail_ratio_quarter;
    records(etaIndex).pg_residual = result.diagnostics.pg_residual;
    records(etaIndex).kkt_residual = result.diagnostics.kkt_residual;
    records(etaIndex).kkt_before_polish = result.diagnostics.kkt_before_polish;
    records(etaIndex).kkt_after_polish = result.diagnostics.kkt_after_polish;
    records(etaIndex).iterations = result.diagnostics.total_iterations;
    records(etaIndex).elapsed_time = result.diagnostics.total_elapsed_time;
    records(etaIndex).tail_mass = result.diagnostics.tail_mass;
    records(etaIndex).tail_max = result.diagnostics.tail_max;
    records(etaIndex).boundary_value_mismatch = ...
        result.diagnostics.boundary_value_mismatch;
    records(etaIndex).polish_status = result.diagnostics.polish_status;
end

fprintf('\neta          Ephys             dEphys      Eaug              minrho     zero_frac   PG/KKT     FFT_tail    iter    time\n');
for j = 1:numberOfEta
    stationarity = max(records(j).pg_residual, records(j).kkt_residual);
    fprintf('%-10.1e  %.12e  %+9.2e  %.12e  %9.2e  %9.2e  %9.2e  %9.2e  %6d  %7.2f\n', ...
        records(j).eta, records(j).physical_energy, ...
        records(j).physical_energy_bias, records(j).augmented_energy, ...
        records(j).min_density, records(j).zero_fraction, stationarity, ...
        records(j).fft_tail_ratio, records(j).iterations, records(j).elapsed_time);
end

biasTolerance = 1e-10 * max(1, abs(Ezero));
baseline_consistent = all([records.physical_energy_bias] >= -biasTolerance);
if ~baseline_consistent
    warning(['Physical energy bias is significantly negative. Do not ' ...
        'interpret the active-set sweep before fixing baseline/solver consistency.']);
end
positiveRecords = records(2:end);
activeReduced = any([positiveRecords.zero_fraction] ...
    < 0.5 * records(1).zero_fraction);
tailImproved = any([positiveRecords.fft_tail_ratio] ...
    < 0.5 * records(1).fft_tail_ratio);
residualCertified = all([positiveRecords.pg_residual] <= 10 * pg_tol);
boxControlled = all(max([positiveRecords.tail_mass], ...
    [positiveRecords.tail_max]) <= 1e-8);
etaZeroHasActiveRegion = records(1).zero_fraction >= 0.05;
if baseline_consistent && etaZeroHasActiveRegion && activeReduced ...
        && tailImproved && residualCertified && boxControlled
    active_set_hypothesis = 'supported';
elseif baseline_consistent && etaZeroHasActiveRegion && activeReduced ...
        && ~tailImproved && residualCertified && boxControlled
    active_set_hypothesis = 'contradicted';
else
    active_set_hypothesis = 'inconclusive';
end
fprintf('\neta=0 disabled/enabled bitwise state/energy: %d / %d\n', ...
    eta_zero_bitwise_state, eta_zero_bitwise_energy);
fprintf('eta=0 continuation consistency: L2 %.3e, Linf %.3e, dE %.3e\n', ...
    eta_zero_continuation_L2, eta_zero_continuation_Linf, ...
    eta_zero_continuation_energy_difference);
fprintf('active-set hypothesis: %s\n', active_set_hypothesis);
fprintf('Evidence: near-zero fraction %.3e -> best %.3e; FFT tail %.3e -> best %.3e.\n', ...
    records(1).zero_fraction, min([positiveRecords.zero_fraction]), ...
    records(1).fft_tail_ratio, min([positiveRecords.fft_tail_ratio]));
fprintf('Evidence: physical energy bias range [%+.3e, %+.3e].\n', ...
    min([positiveRecords.physical_energy_bias]), ...
    max([positiveRecords.physical_energy_bias]));
if strcmp(active_set_hypothesis, 'inconclusive')
    if any(~isfinite([positiveRecords.kkt_residual])) ...
            || ~residualCertified
        likelyRemainingFloor = 'optimizer_residual';
    elseif ~boxControlled
        likelyRemainingFloor = 'box_truncation';
    else
        likelyRemainingFloor = 'nonlinear_aliasing';
    end
    fprintf('Likely remaining floor to check first: %s.\n', likelyRemainingFloor);
else
    likelyRemainingFloor = '';
end

parameters = config.parameters;
parameters.eta_list = eta_list;
archive.parameters = parameters;
archive.regularization = solutions{1}.regularization;
archive.grid = solutions{1}.grid;
archive.solver = config.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = [records.physical_energy].';
archive.physical_energy = [records.physical_energy].';
archive.entropy_value = [records.entropy_value].';
archive.augmented_energy = [records.augmented_energy].';
archive.diagnostics.records = records;
archive.diagnostics.eta_zero_bitwise_state = eta_zero_bitwise_state;
archive.diagnostics.eta_zero_bitwise_energy = eta_zero_bitwise_energy;
archive.diagnostics.eta_zero_continuation_L2 = eta_zero_continuation_L2;
archive.diagnostics.eta_zero_continuation_Linf = eta_zero_continuation_Linf;
archive.diagnostics.eta_zero_continuation_energy_difference = ...
    eta_zero_continuation_energy_difference;
archive.diagnostics.baseline_consistent = baseline_consistent;
archive.diagnostics.active_set_hypothesis = active_set_hypothesis;
archive.diagnostics.likely_remaining_floor = likelyRemainingFloor;
archive.diagnostics.fourier = cellfun( ...
    @(s) s.diagnostics.fourier_tail, solutions, 'UniformOutput', false);
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);

if save_result
    result_file = experiments.SaveResult(archive, 'entropy_sweep', ...
        sprintf('entropy_active_set_N%d_eps_%s.mat', N, ...
        strrep(sprintf('%.0e', epsilon), '-', 'm')), overwrite_existing);
end
if show_plot
    makePlots(solutions, records, eta_list, root, save_result);
end

function record = emptyRecord()
record.eta = NaN;
record.physical_energy = NaN;
record.entropy_value = NaN;
record.augmented_energy = NaN;
record.physical_energy_bias = NaN;
record.density_L2_to_eta0 = NaN;
record.density_Linf_to_eta0 = NaN;
record.min_density = NaN;
record.min_positive_density = NaN;
record.exact_zero_count = NaN;
record.zero_fraction = NaN;
record.zero_tol_list = [];
record.zero_count_by_tol = [];
record.zero_fraction_by_tol = [];
record.support_left_x = NaN;
record.support_right_x = NaN;
record.support_width = NaN;
record.fft_tail_ratio = NaN;
record.pg_residual = NaN;
record.kkt_residual = NaN;
record.kkt_before_polish = NaN;
record.kkt_after_polish = NaN;
record.iterations = NaN;
record.elapsed_time = NaN;
record.tail_mass = NaN;
record.tail_max = NaN;
record.boundary_value_mismatch = NaN;
record.polish_status = '';
end

function makePlots(solutions, records, etaList, root, saveFigures)
labels = arrayfun(@(x) sprintf('eta=%.0e', x), etaList, ...
    'UniformOutput', false);
colors = lines(numel(etaList));
figureFolder = fullfile(root, 'results', 'entropy_sweep');
if saveFigures && ~isfolder(figureFolder)
    mkdir(figureFolder);
end

f1 = figure('Name', 'Entropy density and support endpoint');
tiledlayout(1, 2);
nexttile; hold on;
for j = 1:numel(etaList)
    plot(solutions{j}.grid.x, solutions{j}.rho, 'Color', colors(j, :));
end
xlabel('x'); ylabel('\rho'); grid on; legend(labels, 'Location', 'best');
nexttile; hold on;
base = solutions{1}.diagnostics.active_set;
span = max(2, 0.15 * max(1, base.support_width));
for j = 1:numel(etaList)
    plot(solutions{j}.grid.x, solutions{j}.rho, 'Color', colors(j, :));
end
xlim([base.support_left_x - span, base.support_left_x + span]);
xlabel('x'); ylabel('\rho'); grid on; title('left support endpoint');
saveFigure(f1, figureFolder, '01_density_endpoint.png', saveFigures);

f2 = figure('Name', 'Fourier coefficient decay'); hold on;
for j = 1:numel(etaList)
    fourier = solutions{j}.diagnostics.fourier_tail;
    semilogy(abs(fourier.modes), fourier.rho_hat_abs, '.', ...
        'Color', colors(j, :));
end
xlabel('|k|'); ylabel('abs(rho hat)'); grid on; legend(labels, 'Location', 'best');
set(gca, 'YScale', 'log');
saveFigure(f2, figureFolder, '02_fourier_decay.png', saveFigures);

f3 = figure('Name', 'Active fraction and minimum density');
indices = 1:numel(etaList);
yyaxis left; plot(indices, [records.zero_fraction], 'o-');
ylabel('fraction \rho \leq 10^{-12}');
yyaxis right; plot(indices, [records.min_positive_density], 's--'); hold on;
plot(indices, [records.min_density], 'x-');
set(gca, 'YScale', 'log');
ylabel('min \rho (x), min positive \rho (square)');
xticks(indices); xticklabels(labels); xtickangle(30); grid on;
saveFigure(f3, figureFolder, '03_active_fraction_min_density.png', saveFigures);

f4 = figure('Name', 'Entropy bias and state errors');
semilogy(indices, abs([records.physical_energy_bias]), 'o-', ...
    indices, [records.density_L2_to_eta0], 's-', ...
    indices, [records.density_Linf_to_eta0], '^-');
xticks(indices); xticklabels(labels); xtickangle(30); grid on;
xlim([1, numel(indices)]);
ylabel('absolute diagnostic'); legend('|\Delta E_{phys}|', 'L2 state', ...
    'Linf state', 'Location', 'best');
saveFigure(f4, figureFolder, '04_energy_bias_state_error.png', saveFigures);
end

function saveFigure(handle, folder, name, enabled)
if enabled
    exportgraphics(handle, fullfile(folder, name), 'Resolution', 180);
end
end
