%RUN_PLOT_FOURIER_SPATIAL_CONVERGENCE Paper Fourier mesh-convergence figure.
%
% Fixed (epsilon,sigma) mesh refinement for
%   (i) a harmonic core with boundary-only C-infinity continuation, and
%  (ii) the direct harmonic potential on the periodic box.
% The production solver and the primary spectral comparison are reused
% without modification.
clearvars; clc;

% ======================== fixed experiment ============================
epsilon = 1e-2;
sigma = 1e-2;
beta = 10;
delta = 10;
mass = 1;
L = 32;

% Use the same seven reported meshes for every curve in the spatial-error
% figure.  The compatible checkpoint already contains the two added states.
N_list = [32, 64, 128, 256, 512, 1024, 2048];
N_ref = 8192;
N_list_harm = [32, 64, 128, 256, 512, 1024, 2048];
N_ref_harm = 8192;
N_plot = [64, 128, 256, 512];

R0_fraction = 0.75;
R1_fraction = 0.90;
R0 = R0_fraction * L;
R1 = R1_fraction * L;

s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho)=rho+epsilon';

p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
potential_regularization_label = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma';

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
rate_floor = 1e-13;
reuse_checkpoint = true;
% =====================================================================

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
output_folder = fullfile(root, 'results', 'figures');
if ~isfolder(output_folder)
    mkdir(output_folder);
end
checkpoint_file = fullfile(output_folder, ...
    'mesh_convergence_1d_checkpoint.mat');
spatial_figure_file = fullfile(output_folder, ...
    'mesh_convergence_1d_spatial_error.fig');
spatial_eps_file = fullfile(output_folder, ...
    'mesh_convergence_1d_spatial_error.eps');
fourier_figure_file = fullfile(output_folder, ...
    'mesh_convergence_1d_fourier_decay.fig');
fourier_eps_file = fullfile(output_folder, ...
    'mesh_convergence_1d_fourier_decay.eps');
data_file = fullfile(output_folder, 'mesh_convergence_1d_data.mat');

validateMeshes(N_list, N_ref);
validateMeshes(N_list_harm, N_ref_harm);

template.epsilon = epsilon;
template.sigma = sigma;
template.beta = beta;
template.delta = delta;
template.mass = mass;
template.L = L;
template.s_epsilon = s_epsilon;
template.ds_epsilon = ds_epsilon;
template.d2s_epsilon = d2s_epsilon;
template.fisher_label = fisher_label;
template.p_sigma = p_sigma;
template.dp_sigma = dp_sigma;
template.d2p_sigma = d2p_sigma;
template.potential_regularization_label = ...
    potential_regularization_label;
template.main_pg_switch_tol = main_pg_switch_tol;
template.final_pg_tol = final_pg_tol;

continued_spec.key = 'harmonic_cinf_periodic';
continued_spec.label = ...
    'harmonic core + far-field C-infinity periodic continuation';
continued_spec.V = @(x) model.HarmonicCInfPeriodicPotential( ...
    x, L, R0, R1);
continued_spec.boundary_periodicized = true;
continued_spec.R0 = R0;
continued_spec.R1 = R1;

harmonic_spec.key = 'harmonic';
harmonic_spec.label = 'direct harmonic V(x)=x^2/2';
harmonic_spec.V = @(x) 0.5 * x .^ 2;
harmonic_spec.boundary_periodicized = false;
harmonic_spec.R0 = R0;
harmonic_spec.R1 = R1;

checkpoint = initializeCheckpoint(checkpoint_file, reuse_checkpoint, ...
    template, N_ref, N_ref_harm);

% Solve every intervening factor-two grid so that the N=8192 references
% are reached by Fourier continuation rather than by a single large jump.
continued_solve_N = unique([N_list, powersOfTwoBetween(1024, N_ref), N_ref]);
harmonic_solve_N = unique([N_list_harm, ...
    powersOfTwoBetween(4096, N_ref_harm), N_ref_harm]);

[continued_solutions, checkpoint.continued] = solvePotentialCase( ...
    continued_solve_N, continued_spec, template, checkpoint.continued, ...
    checkpoint_file, checkpoint, 'continued');
[harmonic_solutions, checkpoint.harmonic] = solvePotentialCase( ...
    harmonic_solve_N, harmonic_spec, template, checkpoint.harmonic, ...
    checkpoint_file, checkpoint, 'harmonic');
save(checkpoint_file, 'checkpoint', '-v7.3');

continued_reference = solutionAt(continued_solutions, ...
    continued_solve_N, N_ref);
harmonic_reference = solutionAt(harmonic_solutions, ...
    harmonic_solve_N, N_ref_harm);

[e_res_2, e_res_inf, e_tail_2, e_rho_2, e_E] = ...
    computeErrors(continued_solutions, continued_solve_N, N_list, ...
    continued_reference);
[~, ~, ~, e_rho_2_harm, ~] = computeErrors( ...
    harmonic_solutions, harmonic_solve_N, N_list_harm, ...
    harmonic_reference);

rate_rho2 = rateSeries(e_rho_2, rate_floor);
rate_inf = rateSeries(e_res_inf, rate_floor);
rate_E = rateSeries(e_E, rate_floor);

coeff_decay_data = repmat(struct( ...
    'N', [], 'mode_index', [], 'coefficient_abs', []), ...
    numel(N_plot), 1);
for index = 1:numel(N_plot)
    state = solutionAt(continued_solutions, continued_solve_N, N_plot(index));
    [coefficients, modes] = ...
        src.discretization.ps.FourierCoefficients(state.rho);
    keep = modes >= 0;
    coeff_decay_data(index).N = N_plot(index);
    coeff_decay_data(index).mode_index = modes(keep);
    coeff_decay_data(index).coefficient_abs = abs(coefficients(keep));
end

trapping_potential_labels.continued = continued_spec.label;
trapping_potential_labels.direct_harmonic = harmonic_spec.label;
continued_final_pg = cellfun(@(state) ...
    state.diagnostics.final_pg_residual, continued_solutions);
harmonic_final_pg = cellfun(@(state) ...
    state.diagnostics.final_pg_residual, harmonic_solutions);

save(data_file, 'N_list', 'N_ref', 'e_res_2', 'e_res_inf', ...
    'e_tail_2', 'e_rho_2', 'e_E', 'rate_rho2', 'rate_inf', ...
    'rate_E', 'N_list_harm', 'N_ref_harm', 'e_rho_2_harm', ...
    'coeff_decay_data', 'epsilon', 'sigma', 'beta', 'delta', 'L', ...
    'R0_fraction', 'R1_fraction', 'R0', 'R1', ...
    'trapping_potential_labels', 'continued_final_pg', ...
    'harmonic_final_pg', '-v7.3');

makePaperFigures(N_list, e_rho_2, e_res_inf, e_E, ...
    N_list_harm, e_rho_2_harm, coeff_decay_data, ...
    spatial_figure_file, spatial_eps_file, ...
    fourier_figure_file, fourier_eps_file);
printSummary(N_list, e_rho_2, rate_rho2, e_res_inf, rate_inf, ...
    e_E, rate_E, N_list_harm, e_rho_2_harm, ...
    spatial_figure_file, spatial_eps_file, ...
    fourier_figure_file, fourier_eps_file, data_file, rate_floor);

function [solutions, cache] = solvePotentialCase( ...
    Nvalues, potentialSpec, template, cache, checkpointFile, ...
    wholeCheckpoint, cacheField)
if isempty(cache) || ~isstruct(cache) || ~isfield(cache, 'N_values') ...
        || ~isequal(cache.N_values, Nvalues) ...
        || ~isfield(cache, 'solutions') ...
        || numel(cache.solutions) ~= numel(Nvalues)
    cache.N_values = Nvalues;
    cache.solutions = cell(numel(Nvalues), 1);
end
solutions = cache.solutions;
previous = [];
for index = 1:numel(Nvalues)
    N = Nvalues(index);
    if ~isempty(solutions{index})
        previous = solutions{index};
        fprintf('%s N=%d: reused checkpoint (PG %.3e).\n', ...
            potentialSpec.key, N, ...
            previous.diagnostics.final_pg_residual);
        continue;
    end
    config = makeConfig(N, potentialSpec, template);
    grid = model.SetupGrid1D(config.parameters);
    if isempty(previous)
        rho0 = exp(-grid.x .^ 2);
        rho0 = rho0 * (template.mass ...
            / src.constraints.Mass(rho0, grid.h));
        rho0 = src.constraints.ProjectPositiveConservative( ...
            rho0, template.mass, grid.h, config.solver.projection_tol);
    else
        transferProblem = previous.problem;
        transferProblem.grid = grid;
        transferProblem.plan = src.discretization.ps.Plan1D(grid);
        transferProblem.V = potentialSpec.V(grid.x);
        transferSolver = config.solver;
        transferSolver.splitting = 'legacy_full_gradient';
        rho0 = experiments.TransferState1D(previous.rho, grid, ...
            template.mass, transferProblem, transferSolver);
    end
    fprintf('%s N=%d: solving ...\n', potentialSpec.key, N);
    timer = tic;
    result = experiments.SolveGroundState(config, rho0);
    result.solve_wall_time = toc(timer);
    fprintf('  done %.2f s, E %.15e, PG %.3e\n', ...
        result.solve_wall_time, result.target_energy, ...
        result.diagnostics.final_pg_residual);
    solutions{index} = result;
    previous = result;
    cache.solutions = solutions;
    wholeCheckpoint.(cacheField) = cache;
    checkpoint = wholeCheckpoint; %#ok<NASGU>
    save(checkpointFile, 'checkpoint', '-v7.3');
end
cache.solutions = solutions;
end

function config = makeConfig(N, potentialSpec, template)
config = experiments.DefaultConfig();
config.parameters.beta = template.beta;
config.parameters.delta = template.delta;
config.parameters.mass = template.mass;
config.parameters.L = template.L;
config.parameters.N = N;
config.parameters.epsilon = template.epsilon;

config.fisher_regularization.epsilon = template.epsilon;
config.fisher_regularization.s_epsilon = template.s_epsilon;
config.fisher_regularization.ds_epsilon = template.ds_epsilon;
config.fisher_regularization.d2s_epsilon = template.d2s_epsilon;
config.fisher_regularization.label = template.fisher_label;

config.potential_regularization.name = 'inline_fixed_sigma';
config.potential_regularization.sigma = template.sigma;
config.potential_regularization.p_sigma = template.p_sigma;
config.potential_regularization.dp_sigma = template.dp_sigma;
config.potential_regularization.d2p_sigma = template.d2p_sigma;
config.potential_regularization.label = ...
    template.potential_regularization_label;
config.potential_regularization.prox_type = 'generic_convex';
config.potential_regularization.convexity_tol = 1e-14;

config.trapping_potential.V = potentialSpec.V;
config.trapping_potential.label = potentialSpec.label;
config.trapping_potential.mode = potentialSpec.key;
config.trapping_potential.boundary_periodicized = ...
    potentialSpec.boundary_periodicized;
if potentialSpec.boundary_periodicized
    config.trapping_potential.modification_start = potentialSpec.R0;
    config.trapping_potential.transition_end = potentialSpec.R1;
    config.trapping_potential.reference_V = @(x) 0.5 * x .^ 2;
end

config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.projection_tol = 1e-13;
config.solver.pg_tol = template.main_pg_switch_tol;
config.solver.final_pg_tol = template.final_pg_tol;
config.solver.certification_tol = template.final_pg_tol;
config.solver.residual_check_interval = 10;
config.solver.potential_prox.mass_tol = 1e-12;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.potential_prox.inner_max_iter = 60;
config.solver.max_iter = 200000;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = template.main_pg_switch_tol;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = template.final_pg_tol;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.display = false;
config.output.show_plot = false;
config.output.save_result = false;
end

function [resolvedL2, resolvedLinf, tailL2, totalL2, energyError] = ...
    computeErrors(solutions, solveN, reportN, reference)
resolvedL2 = zeros(size(reportN));
resolvedLinf = zeros(size(reportN));
tailL2 = zeros(size(reportN));
totalL2 = zeros(size(reportN));
energyError = zeros(size(reportN));
for index = 1:numel(reportN)
    state = solutionAt(solutions, solveN, reportN(index));
    comparison = src.diagnostics.SpectralStateComparison( ...
        state.rho, state.grid, reference.rho, reference.grid);
    resolvedL2(index) = comparison.resolved_L2_error;
    resolvedLinf(index) = comparison.resolved_Linf_error;
    tailL2(index) = comparison.reference_tail_L2;
    totalL2(index) = comparison.total_spectral_L2_error;
    energyError(index) = abs(state.target_energy - reference.target_energy);
end
end

function state = solutionAt(solutions, Nvalues, N)
index = find(Nvalues == N, 1);
if isempty(index) || isempty(solutions{index})
    error('Missing computed state at N=%d.', N);
end
state = solutions{index};
end

function rates = rateSeries(errors, floorValue)
rates = NaN(size(errors));
for index = 2:numel(errors)
    if isfinite(errors(index - 1)) && isfinite(errors(index)) ...
            && errors(index - 1) > floorValue ...
            && errors(index) > floorValue && errors(index) > 0
        rates(index) = log2(errors(index - 1) / errors(index));
    end
end
end

function makePaperFigures(Nvalues, stateL2, stateLinf, energyError, ...
    Nharm, stateL2Harm, coefficientData, spatialFigFile, ...
    spatialEpsFile, fourierFigFile, fourierEpsFile)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
spatialFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80, 80, 980, 780], ...
    'Name', 'Fourier spatial convergence');
axisA = axes(spatialFigure);
loglog(axisA, Nvalues, stateL2, '-o', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'continued harmonic: $e_{\rho,2}$');
hold(axisA, 'on');
loglog(axisA, Nvalues, stateLinf, '--s', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'continued harmonic: $e_{\rho,\infty}$');
loglog(axisA, Nvalues, energyError, '-.^', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'continued harmonic: $e_E$');
loglog(axisA, Nharm, stateL2Harm, ':d', ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'direct harmonic: $e_{\rho,2}$');
% The high-resolution direct-harmonic tail is second order.  Plot a
% parallel O(N^{-2}) guide over its final three mesh levels, with a
% small vertical offset so that the reference remains visible.
tailStart = max(1, numel(Nharm) - 2);
referenceN = Nharm(tailStart:end);
referenceError = 2.5 * stateL2Harm(tailStart) ...
    .* (referenceN ./ referenceN(1)).^(-2);
loglog(axisA, referenceN, referenceError, 'k--', ...
    'LineWidth', lineWidth, ...
    'DisplayName', '$\mathcal{O}(N^{-2})$');
grid(axisA, 'on'); box(axisA, 'on');
xlabel(axisA, '$N$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(axisA, 'error', 'Interpreter', 'latex', 'FontSize', fontSize);
title(axisA, '(a) Fourier spatial convergence', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(axisA, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(axisA, 'FontSize', fontSize, 'LineWidth', 1);
% Save an explicitly visible FIG. Otherwise openfig/double-click can load
% the object successfully while leaving its window hidden.
set(spatialFigure, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(spatialFigure, spatialFigFile);
set(spatialFigure, 'Visible', 'off');
print(spatialFigure, spatialEpsFile, '-depsc2', '-painters');
close(spatialFigure);

fourierFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80, 80, 980, 780], ...
    'Name', 'Fourier coefficient decay');
axisB = axes(fourierFigure);
lineStyles = {'-o', '--s', '-.^', ':d'};
set(axisB, 'YScale', 'log');
hold(axisB, 'on');
for index = 1:numel(coefficientData)
    modes = coefficientData(index).mode_index;
    amplitudes = coefficientData(index).coefficient_abs;
    markerIndices = unique(round(linspace(1, numel(modes), ...
        min(24, numel(modes)))));
    semilogy(axisB, modes, max(amplitudes, realmin), ...
        lineStyles{index}, 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, 'MarkerIndices', markerIndices, ...
        'DisplayName', sprintf('$N=%d$', coefficientData(index).N));
end
grid(axisB, 'on'); box(axisB, 'on');
xlabel(axisB, 'mode index $|k|$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
ylabel(axisB, '$|\widehat{\rho}_k|$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
title(axisB, '(b) Fourier coefficient decay', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(axisB, 'Location', 'northeast', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(axisB, 'FontSize', fontSize, 'LineWidth', 1);
set(fourierFigure, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(fourierFigure, fourierFigFile);
set(fourierFigure, 'Visible', 'off');
print(fourierFigure, fourierEpsFile, '-depsc2', '-painters');
close(fourierFigure);
end

function printSummary(Nvalues, stateL2, rateL2, stateLinf, rateInf, ...
    energyError, rateEnergy, Nharm, stateL2Harm, ...
    spatialFigFile, spatialEpsFile, fourierFigFile, fourierEpsFile, ...
    dataFile, floorValue)
fprintf('\n============================================================\n');
fprintf('Continued harmonic Fourier spatial convergence\n');
fprintf('============================================================\n');
fprintf(' N       e_rho_2    rate      e_res_inf  rate      e_E        rate\n');
for index = 1:numel(Nvalues)
    fprintf('%4d   %.3e  %7s   %.3e  %7s   %.3e  %7s\n', ...
        Nvalues(index), stateL2(index), ...
        rateLabel(rateL2(index), stateL2(index), floorValue), ...
        stateLinf(index), ...
        rateLabel(rateInf(index), stateLinf(index), floorValue), ...
        energyError(index), ...
        rateLabel(rateEnergy(index), energyError(index), floorValue));
end
fprintf('\nDirect harmonic comparison\n');
fprintf(' N       e_rho_2_harm\n');
for index = 1:numel(Nharm)
    fprintf('%4d   %.3e\n', Nharm(index), stateL2Harm(index));
end
fprintf('\nSaved outputs\n');
fprintf('  %s\n', spatialFigFile);
fprintf('  %s\n', spatialEpsFile);
fprintf('  %s\n', fourierFigFile);
fprintf('  %s\n', fourierEpsFile);
fprintf('  %s\n', dataFile);
end

function label = rateLabel(rate, errorValue, floorValue)
if errorValue <= floorValue
    label = 'floor';
elseif isfinite(rate)
    label = sprintf('%.2f', rate);
else
    label = '--';
end
end

function validateMeshes(Nvalues, Nreference)
if any(mod(Nvalues, 2) ~= 0) || mod(Nreference, 2) ~= 0 ...
        || Nreference <= max(Nvalues) ...
        || any(mod(Nreference, Nvalues) ~= 0) ...
        || any(diff(Nvalues) <= 0)
    error('Use strictly increasing nested even Fourier grids.');
end
end

function values = powersOfTwoBetween(firstValue, lastValue)
if firstValue > lastValue
    values = [];
    return;
end
values = firstValue;
while values(end) < lastValue
    values(end + 1) = 2 * values(end); %#ok<AGROW>
end
values = values(values <= lastValue);
end

function checkpoint = initializeCheckpoint(filename, reuse, ...
    template, Nreference, NreferenceHarmonic)
checkpoint.signature = [template.epsilon, template.sigma, ...
    template.beta, template.delta, template.mass, template.L, ...
    Nreference, NreferenceHarmonic];
checkpoint.continued = [];
checkpoint.harmonic = [];
if ~reuse || ~isfile(filename)
    return;
end
loaded = load(filename, 'checkpoint');
if isfield(loaded, 'checkpoint') ...
        && isfield(loaded.checkpoint, 'signature') ...
        && isequal(loaded.checkpoint.signature, checkpoint.signature)
    checkpoint = loaded.checkpoint;
    fprintf('Reusing compatible Fourier-convergence checkpoint.\n');
else
    warning('Ignoring incompatible Fourier-convergence checkpoint.');
end
end
