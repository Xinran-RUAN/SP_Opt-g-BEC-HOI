%RUN_POTENTIAL_REGULARIZATION_EFFECTS_L32 Final potential-smoothing study.
%
% This experiment recomputes an independent unsmoothed baseline and four
% fixed-sigma minimizers on the same L=32 Fourier grid. No L=8 state or
% checkpoint is eligible for reuse.
clearvars; clc;

% ======================== final experiment ============================
epsilon = 1e-2;
L = 32;
N = 4096;
beta = 10;
delta = 10;
mass = 1;
power_list = 1:4;
sigma_list = epsilon .^ power_list;

s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho)=rho+epsilon';
potential_definition = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma';
trapping_potential_definition = 'V(x)=x^2/2';

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-10;
reuse_checkpoint = true;
% =====================================================================

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
output_folder = fullfile(root, 'results', 'potential_sigma');
if ~isfolder(output_folder)
    mkdir(output_folder);
end
checkpoint_file = fullfile(output_folder, ...
    'regularization_sigma_L32_checkpoint.mat');
bias_figure_file = fullfile(output_folder, ...
    'regularization_sigma_L32_bias.fig');
bias_eps_file = fullfile(output_folder, ...
    'regularization_sigma_L32_bias.eps');
transition_figure_file = fullfile(output_folder, ...
    'regularization_sigma_L32_transition.fig');
transition_eps_file = fullfile(output_folder, ...
    'regularization_sigma_L32_transition.eps');
data_file = fullfile(output_folder, 'regularization_sigma_L32_data.mat');

signature.epsilon = epsilon;
signature.L = L;
signature.N = N;
signature.beta = beta;
signature.delta = delta;
signature.mass = mass;
signature.power_list = power_list;
signature.sigma_list = sigma_list;
signature.fisher_definition = fisher_label;
signature.potential_definition = potential_definition;
signature.trapping_potential_definition = trapping_potential_definition;
checkpoint = loadCompatibleCheckpoint( ...
    checkpoint_file, signature, reuse_checkpoint);

base_config = makeBaseConfig(epsilon, L, N, beta, delta, mass, ...
    s_epsilon, ds_epsilon, d2s_epsilon, fisher_label, ...
    main_pg_switch_tol, final_pg_tol);
grid = model.SetupGrid1D(base_config.parameters);
rho_initial = exp(-grid.x .^ 2) / sqrt(pi);
rho_initial = rho_initial * (mass ...
    / src.constraints.Mass(rho_initial, grid.h));
rho_initial = src.constraints.ProjectPositiveConservative( ...
    rho_initial, mass, grid.h, base_config.solver.projection_tol);

if isempty(checkpoint.baseline)
    baseline_config = setLinearPotential(base_config);
    fprintf('Solving independent L=32, N=%d unsmoothed baseline ...\n', N);
    checkpoint.baseline = experiments.SolveGroundState( ...
        baseline_config, rho_initial);
    saveCheckpoint(checkpoint_file, checkpoint);
else
    fprintf('Reused compatible L=32, N=%d baseline checkpoint.\n', N);
end
baseline_result = checkpoint.baseline;

for index = 1:numel(sigma_list)
    sigma = sigma_list(index);
    if ~isempty(checkpoint.solutions{index})
        fprintf('Reused compatible sigma=%.0e checkpoint.\n', sigma);
        continue;
    end
    variant_config = setSmoothedPotential(base_config, sigma, ...
        potential_definition);
    fprintf('Solving L=32, N=%d, sigma=%.0e ...\n', N, sigma);
    % Sigma continuation is used only as a warm start. Every state is
    % subsequently minimized and certified for its own convex objective.
    rho_start = rho_initial;
    if index > 1 && ~isempty(checkpoint.solutions{index - 1})
        rho_start = checkpoint.solutions{index - 1}.rho;
    end
    checkpoint.solutions{index} = experiments.SolveGroundState( ...
        variant_config, rho_start);
    saveCheckpoint(checkpoint_file, checkpoint);
end
solutions = checkpoint.solutions;

records = repmat(src.diagnostics.PotentialStudyDiagnostics( ...
    solutions{1}, baseline_result), numel(sigma_list), 1);
for index = 2:numel(sigma_list)
    records(index) = src.diagnostics.PotentialStudyDiagnostics( ...
        solutions{index}, baseline_result);
end
for index = 1:numel(records)
    records(index).min_density_over_sigma = ...
        records(index).min_density / records(index).sigma;
end

% Use the code's diagnostic definition verbatim. For L=32 this gives
% ell_tail=min(2,0.1L)=2 and D_L^tail={|x|>=30}.
tail_definition.ell_tail = min(2, 0.1 * L);
tail_definition.threshold = L - tail_definition.ell_tail;
tail_definition.expression = ...
    'D_L^tail={x in [-L,L): |x| >= L-min(2,0.1L)}';
tail_definition.discrete_rule = ...
    'h*sum(rho(abs(x)>=L-min(2,0.1L)))';

baseline_metadata.target_energy = baseline_result.target_energy;
baseline_metadata.baseline_energy = baseline_result.baseline_energy;
baseline_metadata.final_pg = baseline_result.diagnostics.final_pg_residual;
baseline_metadata.exact_zero_count = ...
    baseline_result.diagnostics.exact_zero_count;

q_sigma_profiles = cell(numel(sigma_list), 1);
for index = 1:numel(sigma_list)
    sigma = sigma_list(index);
    q_sigma_profiles{index} = solutions{index}.rho ...
        ./ hypot(solutions{index}.rho, sigma);
end

printTables(records, baseline_metadata, epsilon, L, N, beta, delta, ...
    tail_definition);
makeFigures(grid.x, records, q_sigma_profiles, ...
    bias_figure_file, bias_eps_file, ...
    transition_figure_file, transition_eps_file);

rho_baseline = baseline_result.rho;
rho_sigma = cellfun(@(state) state.rho, solutions, ...
    'UniformOutput', false);
target_energy = [records.target_energy].';
baseline_energy = [records.baseline_energy].';
B_E = [records.baseline_energy_bias].';
B_rho = [records.density_L2_bias].';
m_tail = [records.tail_mass].';
min_rho = [records.min_density].';
final_PG = [records.final_pg_residual].';
exact_zero_count = [records.exact_zero_count].';
m_tail_over_sigma = [records.tail_mass_over_sigma].';
min_rho_over_sigma = [records.min_density_over_sigma].';
x = grid.x;
h = grid.h;
metadata.epsilon = epsilon;
metadata.L = L;
metadata.N = N;
metadata.h = h;
metadata.beta = beta;
metadata.delta = delta;
metadata.mass = mass;
metadata.power_list = power_list;
metadata.sigma_list = sigma_list;
metadata.fisher_definition = fisher_label;
metadata.potential_definition = potential_definition;
metadata.trapping_potential_definition = trapping_potential_definition;
metadata.tail_definition = tail_definition;
metadata.checkpoint_signature = signature;
save(data_file, 'metadata', 'epsilon', 'L', 'N', 'h', 'beta', ...
    'delta', 'mass', 'power_list', 'sigma_list', 'x', ...
    'rho_baseline', 'rho_sigma', ...
    'target_energy', 'baseline_energy', 'B_E', 'B_rho', 'm_tail', ...
    'min_rho', 'final_PG', 'exact_zero_count', ...
    'm_tail_over_sigma', 'min_rho_over_sigma', 'q_sigma_profiles', ...
    'tail_definition', 'baseline_metadata', 'records', '-v7.3');

fprintf('\nSaved outputs\n');
fprintf('  %s\n', bias_figure_file);
fprintf('  %s\n', bias_eps_file);
fprintf('  %s\n', transition_figure_file);
fprintf('  %s\n', transition_eps_file);
fprintf('  %s\n', data_file);

function config = makeBaseConfig(epsilon, L, N, beta, delta, mass, ...
    s, ds, d2s, fisherLabel, mainTolerance, finalTolerance)
config = experiments.DefaultConfig();
config.parameters.epsilon = epsilon;
config.parameters.L = L;
config.parameters.N = N;
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.mass = mass;
config.trapping_potential.V = @(x) 0.5 * x .^ 2;
config.trapping_potential.label = 'V(x)=x^2/2';
config.fisher_regularization.epsilon = epsilon;
config.fisher_regularization.s_epsilon = s;
config.fisher_regularization.ds_epsilon = ds;
config.fisher_regularization.d2s_epsilon = d2s;
config.fisher_regularization.label = fisherLabel;
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.projection_tol = 1e-13;
config.solver.pg_tol = mainTolerance;
config.solver.final_pg_tol = finalTolerance;
config.solver.certification_tol = finalTolerance;
config.solver.residual_check_interval = 10;
config.solver.potential_prox.mass_tol = 1e-12;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.potential_prox.inner_max_iter = 80;
config.solver.max_iter = 200000;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = mainTolerance;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = finalTolerance;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.display = false;
config.output.show_plot = false;
config.output.save_result = false;
end

function config = setLinearPotential(config)
config.potential_regularization.name = 'inline_linear';
config.potential_regularization.sigma = 0;
config.potential_regularization.p_sigma = @(rho) rho;
config.potential_regularization.dp_sigma = @(rho) ones(size(rho));
config.potential_regularization.d2p_sigma = @(rho) zeros(size(rho));
config.potential_regularization.label = 'p_0(rho)=rho';
config.potential_regularization.prox_type = 'linear';
end

function config = setSmoothedPotential(config, sigma, label)
config.potential_regularization.name = 'inline_p_sigma';
config.potential_regularization.sigma = sigma;
config.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
config.potential_regularization.dp_sigma = @(rho) ...
    rho ./ hypot(rho, sigma);
config.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
config.potential_regularization.label = label;
config.potential_regularization.prox_type = 'generic_convex';
config.potential_regularization.convexity_tol = 1e-14;
end

function checkpoint = loadCompatibleCheckpoint(filename, signature, reuse)
checkpoint.signature = signature;
checkpoint.baseline = [];
checkpoint.solutions = cell(numel(signature.sigma_list), 1);
if ~reuse || ~isfile(filename)
    return;
end
loaded = load(filename, 'checkpoint');
if isfield(loaded, 'checkpoint') ...
        && isfield(loaded.checkpoint, 'signature') ...
        && isequaln(loaded.checkpoint.signature, signature) ...
        && isfield(loaded.checkpoint, 'baseline') ...
        && isfield(loaded.checkpoint, 'solutions') ...
        && numel(loaded.checkpoint.solutions) == numel(signature.sigma_list)
    checkpoint = loaded.checkpoint;
else
    warning(['Ignoring incompatible sigma checkpoint. Compatibility ' ...
        'requires epsilon, sigma, L, N, beta, delta, Fisher, potential, ' ...
        'and trapping-potential definitions to agree.']);
end
end

function saveCheckpoint(filename, checkpoint)
save(filename, 'checkpoint', '-v7.3');
end

function printTables(records, baseline, epsilon, L, N, beta, delta, tail)
fprintf('\n============================================================\n');
fprintf('Potential smoothing / Regularization effects\n');
fprintf('epsilon=%.0e, L=%g, N=%d, h=%.6g, beta=%g, delta=%g\n', ...
    epsilon, L, N, 2 * L / N, beta, delta);
fprintf('Baseline E_eps,0 = %.15e, PG = %.3e, exact0 = %d\n', ...
    baseline.baseline_energy, baseline.final_pg, ...
    baseline.exact_zero_count);
fprintf('Tail: ell_tail=%.6g, |x| >= %.6g\n', ...
    tail.ell_tail, tail.threshold);
fprintf('============================================================\n');
fprintf([' sigma      B_E          B_rho        m_tail       min_rho      ' ...
    'PG          exact0\n']);
for index = 1:numel(records)
    r = records(index);
    fprintf(' %.0e   %+.6e  %.6e  %.6e  %.6e  %.3e  %d\n', ...
        r.sigma, r.baseline_energy_bias, r.density_L2_bias, ...
        r.tail_mass, r.min_density, r.final_pg_residual, ...
        r.exact_zero_count);
end
fprintf('\nRatio diagnostics\n');
fprintf(' sigma      m_tail/sigma   min_rho/sigma\n');
for index = 1:numel(records)
    fprintf(' %.0e   %.6e     %.6e\n', records(index).sigma, ...
        records(index).tail_mass_over_sigma, ...
        records(index).min_density_over_sigma);
end
end

function makeFigures(x, records, qProfiles, biasFigFile, biasEpsFile, ...
    transitionFigFile, transitionEpsFile)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
colors = lines(numel(records));
sigmaValues = [records.sigma];

% Figure 1: regularization bias and tail mass.  This is deliberately a
% stand-alone paper figure rather than one panel of a tiled layout.
fBias = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 650], ...
    'Name', 'Potential smoothing effects');
ax1 = axes('Parent', fBias);
loglog(ax1, sigmaValues, [records.baseline_energy_bias], '-o', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', '$B_E$');
hold(ax1, 'on');
loglog(ax1, sigmaValues, [records.density_L2_bias], '--s', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', '$B_\rho$');
loglog(ax1, sigmaValues, [records.min_density], '-.^', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', '$\min\rho$');
referenceScale = sqrt(records(1).baseline_energy_bias ...
    * records(1).density_L2_bias) / sigmaValues(1);
referenceValues = referenceScale .* sigmaValues;
loglog(ax1, sigmaValues, referenceValues, 'k:', ...
    'LineWidth', lineWidth, ...
    'DisplayName', '$\mathcal{O}(\sigma)$');
set(ax1, 'XDir', 'reverse', 'FontSize', fontSize, 'LineWidth', 1);
xMarginFactor = 2;
xlim(ax1, [min(sigmaValues) / xMarginFactor, ...
    max(sigmaValues) * xMarginFactor]);
allBiasValues = [[records.baseline_energy_bias], ...
    [records.density_L2_bias], [records.min_density], referenceValues];
ylim(ax1, [min(allBiasValues) / 2, max(allBiasValues) * 5]);
grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, '$\sigma$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(ax1, 'error', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
title(ax1, 'Potential-smoothing effects', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(ax1, 'Location', 'best', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(fBias, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(fBias, biasFigFile);
set(fBias, 'Visible', 'off', 'Renderer', 'painters');
print(fBias, biasEpsFile, '-depsc2', '-vector');
close(fBias);

% Figure 2: show the whole right half-domain and add a reproducible inset
% around the actual q_sigma transition.  The inset is computed from the
% stored profiles; no solver state is altered or recomputed.
fTransition = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 650], ...
    'Name', 'Regularized transition profiles');
ax2 = axes('Parent', fTransition, 'Position', [0.12, 0.13, 0.82, 0.79]);
lineStyles = {'-', '--', '-.', ':'};
hold(ax2, 'on');
right = x >= 0;
[xLower, xUpper] = transitionWindow(x, qProfiles);
for index = 1:numel(records)
    indices = find(right);
    markerIndices = unique(round(linspace(1, numel(indices), ...
        min(20, numel(indices)))));
    plot(ax2, x(indices), qProfiles{index}(indices), ...
        lineStyles{index}, 'Color', colors(index, :), ...
        'LineWidth', lineWidth, 'Marker', markerFor(index), ...
        'MarkerSize', markerSize, 'MarkerIndices', markerIndices, ...
        'DisplayName', sprintf('$\\sigma=10^{%d}$', ...
        round(log10(records(index).sigma))));
end
set(ax2, 'FontSize', fontSize, 'LineWidth', 1);
grid(ax2, 'on'); box(ax2, 'on');
xlim(ax2, [0, max(x(right))]); ylim(ax2, [-0.02, 1.02]);
xlabel(ax2, '$x$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(ax2, '$q_\sigma(x)$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
title(ax2, 'Regularized transition profiles', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(ax2, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', fontSize);

% Mark and enlarge the transition interval.  The inset repeats the same
% stored data without applying a plotting floor or modifying rho.
rectangle(ax2, 'Position', [xLower, 0, xUpper - xLower, 1], ...
    'EdgeColor', [0.30, 0.30, 0.30], 'LineStyle', ':', ...
    'LineWidth', 1.25, 'HandleVisibility', 'off');
axInset = axes('Parent', fTransition, ...
    'Position', [0.53, 0.46, 0.38, 0.36]);
hold(axInset, 'on');
for index = 1:numel(records)
    indices = find(right & x >= xLower & x <= xUpper);
    markerIndices = unique(round(linspace(1, numel(indices), ...
        min(10, numel(indices)))));
    plot(axInset, x(indices), qProfiles{index}(indices), ...
        lineStyles{index}, 'Color', colors(index, :), ...
        'LineWidth', lineWidth, 'Marker', markerFor(index), ...
        'MarkerSize', 6, 'MarkerIndices', markerIndices, ...
        'HandleVisibility', 'off');
end
set(axInset, 'FontSize', 14, 'LineWidth', 1);
grid(axInset, 'on'); box(axInset, 'on');
xlim(axInset, [xLower, xUpper]); ylim(axInset, [-0.02, 1.02]);
xlabel(axInset, '$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(axInset, '$q_\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
title(axInset, 'transition-layer zoom', ...
    'FontSize', 14, 'FontWeight', 'normal');

set(fTransition, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(fTransition, transitionFigFile);
set(fTransition, 'Visible', 'off', 'Renderer', 'painters');
print(fTransition, transitionEpsFile, '-depsc2', '-vector');
close(fTransition);
end

function [lower, upper] = transitionWindow(x, profiles)
positive = x >= 0;
locations = [];
for index = 1:numel(profiles)
    % Use the principal q_sigma transition rather than the long, tiny-q
    % far-field tail; otherwise the inset would span almost the full box.
    transition = positive & profiles{index} >= 0.10 ...
        & profiles{index} <= 0.90;
    locations = [locations; x(transition)]; %#ok<AGROW>
end
if isempty(locations)
    lower = 0;
    upper = min(max(x), 0.25 * (max(x) - min(x)));
else
    padding = max(0.5, 0.1 * (max(locations) - min(locations)));
    lower = max(0, min(locations) - padding);
    upper = min(max(x), max(locations) + padding);
end
end

function marker = markerFor(index)
markers = {'o', 's', '^', 'd'};
marker = markers{index};
end
