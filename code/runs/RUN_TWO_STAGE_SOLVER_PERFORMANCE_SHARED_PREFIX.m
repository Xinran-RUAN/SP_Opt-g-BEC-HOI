%RUN_TWO_STAGE_SOLVER_PERFORMANCE_SHARED_PREFIX Shared-prefix Fig. 5.5.
%
% One uninterrupted FISTA-CD run supplies both the common prefix and the
% FISTA-continuation branch, preserving its momentum, accepted step, and
% restart state.  The production Newton--PCG refinement starts from the
% captured switch state.  No production solver algorithm is modified.

if ~exist('reuse_certified_branch_data', 'var')
    reuse_certified_branch_data = true;
end
clearvars -except reuse_certified_branch_data;
clc;

% =========================== fixed setup =============================
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 4096;
epsilon = 1e-2;
sigma = 1e-4;
R0_fraction = 0.75;
R1_fraction = 0.90;
fista_plotting_cutoff = 50000;

r_epsilon = @(rho) rho + epsilon;
dr_epsilon = @(rho) ones(size(rho));
d2r_epsilon = @(rho) zeros(size(rho));
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
% =====================================================================

projectRoot = fileparts(fileparts(mfilename('fullpath')));
repositoryRoot = fileparts(projectRoot);
addpath(projectRoot);
startup_HOI();

resultFolder = fullfile(projectRoot, 'results', 'solver_performance');
figureFolder = fullfile(projectRoot, 'figs');
manuscriptFigureFolder = fullfile(repositoryRoot, 'manuscript', 'figs');
if ~isfolder(resultFolder), mkdir(resultFolder); end
if ~isfolder(figureFolder), mkdir(figureFolder); end

dataFile = fullfile(resultFolder, ...
    'two_stage_vs_fista_shared_prefix.mat');
legacyFile = fullfile(resultFolder, 'two_stage_vs_fista_only.mat');
energyBase = fullfile(figureFolder, ...
    'solver_energy_history_shared_prefix');
residualBase = fullfile(figureFolder, ...
    'solver_residual_history_shared_prefix');
diagnosticBase = fullfile(resultFolder, ...
    'shared_prefix_iteration_overlap');

R0 = R0_fraction * L;
R1 = R1_fraction * L;
config = makeConfig(beta, delta, mass, L, N, R0, R1, ...
    epsilon, sigma, r_epsilon, dr_epsilon, d2r_epsilon, ...
    p_sigma, dp_sigma, d2p_sigma);
grid = model.SetupGrid1D(config.parameters);
rho0 = model.InitialDensity(grid, mass, 'semismooth');

source = 'fresh shared-prefix branch run';
if reuse_certified_branch_data && isfile(legacyFile)
    legacy = load(legacyFile);
    validateLegacyArchive(legacy, epsilon, sigma, L, N, grid.h);
    fprintf(['Reusing the certified branch histories.  The independent ' ...
        'two-stage prefix timing is discarded.\n']);
    fistaResult = legacy.fistaOnlyResult;
    newtonHistory = legacy.twoStageResult.history.polish;
    newtonFinalRho = legacy.twoStageResult.rho;
    newtonFinalEnergy = legacy.twoStageResult.target_energy;
    newtonFinalPg = legacy.twoStageResult.diagnostics.final_pg_residual;
    newtonFinalKkt = legacy.twoStageResult.diagnostics.final_kkt_residual;
    newtonIterations = legacy.twoStageResult.diagnostics.polish_iterations;
    newtonElapsedTime = legacy.twoStageResult.diagnostics.polish_elapsed_time;
    newtonPcgZ = legacy.twoStageResult.diagnostics.total_pcg_z_iterations;
    newtonPcgW = legacy.twoStageResult.diagnostics.total_pcg_w_iterations;
    newtonPcgMax = legacy.twoStageResult.diagnostics.max_pcg_iterations;
    canonicalSolver = legacy.twoStageResult.solver;
    source = ['certified legacy branches, reorganized with the sole ' ...
        'FISTA-only prefix'];
else
    fprintf('Running uninterrupted FISTA prefix/continuation branch ...\n');
    branchConfig = config;
    branchConfig.solver.polish_mode = 'none';
    branchConfig.solver.switch.stop_at_energy_handoff = false;
    branchConfig.solver.history.capture_handoff_state = true;
    branchConfig.solver.max_iter = fista_plotting_cutoff;
    fistaResult = experiments.SolveGroundState(branchConfig, rho0);

    fprintf('Running Newton--PCG from the captured switch state ...\n');
    switchState = fistaResult.diagnostics.handoff_rho;
    polishOptions = makePolishOptions(fistaResult.solver);
    newtonResult = src.solvers.PolishKKT( ...
        switchState, fistaResult.problem, polishOptions);
    if newtonResult.failed || ~newtonResult.polish_converged
        warning('Shared-prefix Newton branch was not fully certified: %s', ...
            newtonResult.status);
    end
    newtonHistory = newtonResult.history;
    newtonFinalRho = newtonResult.rho;
    newtonFinalEnergy = newtonResult.energy;
    newtonFinalPg = newtonResult.full_pg_residual;
    newtonFinalKkt = newtonResult.kkt_residual;
    newtonIterations = newtonResult.iterations;
    newtonElapsedTime = newtonResult.elapsed_time;
    newtonPcgZ = newtonResult.total_pcg_z_iterations;
    newtonPcgW = newtonResult.total_pcg_w_iterations;
    newtonPcgMax = newtonResult.max_pcg_iterations;
    canonicalSolver = config.solver;
end

problem = fistaResult.problem;
solver = fistaResult.solver;
mainHistory = fistaResult.history.main;
switchIteration = fistaResult.diagnostics.handoff_iteration;
assert(isfinite(switchIteration) && switchIteration >= 1, ...
    'The common FISTA run did not detect the production switch.');
assert(numel(mainHistory.iteration) >= switchIteration, ...
    'The FISTA continuation history ends before the switch.');
switchState = fistaResult.diagnostics.handoff_rho;
assert(~isempty(switchState), 'The switch state was not archived.');

initialEnergy = src.discretization.ps.Energy(rho0, problem);
initialGradient = src.discretization.ps.Gradient(rho0, problem);
initialPg = src.solvers.FullGradientMapping( ...
    rho0, initialGradient, problem, solver);

% Common FISTA prefix, including iteration zero.
common_prefix.iteration = (0:switchIteration).';
common_prefix.time = [0; mainHistory.elapsed_time(1:switchIteration)];
common_prefix.energy = [initialEnergy; ...
    mainHistory.augmented_energy(1:switchIteration)];
common_prefix.pg = [initialPg; ...
    mainHistory.full_pg_residual(1:switchIteration)];
common_prefix.accepted_L = [NaN; ...
    mainHistory.accepted_L(1:switchIteration)];
common_prefix.backtracks = [NaN; ...
    mainHistory.backtracks(1:switchIteration)];
common_prefix.restart = [false; ...
    mainHistory.restart(1:switchIteration)];

switchTime = common_prefix.time(end);
switchEnergy = common_prefix.energy(end);
switchPg = common_prefix.pg(end);

% Incremental FISTA continuation history, including the switch state at dt=0.
continuationIndex = (switchIteration + 1):numel(mainHistory.iteration);
fista_continuation.dt = [0; ...
    mainHistory.elapsed_time(continuationIndex) - switchTime];
fista_continuation.energy = [switchEnergy; ...
    mainHistory.augmented_energy(continuationIndex)];
fista_continuation.pg = [switchPg; ...
    mainHistory.full_pg_residual(continuationIndex)];
fista_continuation.iteration = [switchIteration; ...
    mainHistory.iteration(continuationIndex)];

% Incremental Newton history, including the switch state at dt=0.
assert(isfield(newtonHistory, 'elapsed_time') ...
    && isfield(newtonHistory, 'energy') ...
    && isfield(newtonHistory, 'full_pg_residual'), ...
    'The Newton branch does not contain the required production history.');
newton_branch.dt = [0; newtonHistory.elapsed_time(:)];
newton_branch.energy = [switchEnergy; newtonHistory.energy(:)];
newton_branch.pg = [switchPg; newtonHistory.full_pg_residual(:)];
newton_branch.iteration = [0; newtonHistory.iteration(:)];
% The production Newton timer includes its final standard KKT packaging,
% whereas the per-iteration history ends at the last accepted iterate.
% Append that final solver timestamp without changing the numerical state.
if newtonElapsedTime > newton_branch.dt(end)
    newton_branch.dt(end + 1, 1) = newtonElapsedTime;
    newton_branch.energy(end + 1, 1) = newtonFinalEnergy;
    newton_branch.pg(end + 1, 1) = newtonFinalPg;
    newton_branch.iteration(end + 1, 1) = newtonIterations;
else
    newton_branch.dt(end) = newtonElapsedTime;
    newton_branch.energy(end) = newtonFinalEnergy;
    newton_branch.pg(end) = newtonFinalPg;
end
newton_branch.final_rho = newtonFinalRho;
newton_branch.final_energy = newtonFinalEnergy;
newton_branch.final_pg = newtonFinalPg;
newton_branch.final_kkt = newtonFinalKkt;
newton_branch.iterations = newtonIterations;
newton_branch.elapsed_time = newtonElapsedTime;

% Formal plotted trajectories.  The branch-initial switch point is not
% duplicated when the incremental histories are appended.
plotted_fista.time = [common_prefix.time; ...
    switchTime + fista_continuation.dt(2:end)];
plotted_fista.energy = [common_prefix.energy; ...
    fista_continuation.energy(2:end)];
plotted_fista.pg = [common_prefix.pg; fista_continuation.pg(2:end)];

plotted_two_stage.time = [common_prefix.time; ...
    switchTime + newton_branch.dt(2:end)];
plotted_two_stage.energy = [common_prefix.energy; ...
    newton_branch.energy(2:end)];
plotted_two_stage.pg = [common_prefix.pg; newton_branch.pg(2:end)];

E_ref = newtonFinalEnergy;
energyPlotFloor = max(1e-16, 10 * eps(max(1, abs(E_ref))));
plotted_fista.energy_error_raw = abs(plotted_fista.energy - E_ref);
plotted_fista.energy_error_plot = max( ...
    plotted_fista.energy_error_raw, energyPlotFloor);
plotted_two_stage.energy_error_raw = abs( ...
    plotted_two_stage.energy - E_ref);
plotted_two_stage.energy_error_plot = max( ...
    plotted_two_stage.energy_error_raw, energyPlotFloor);

prefixLength = numel(common_prefix.time);
checks.max_prefix_time_difference = max(abs( ...
    plotted_fista.time(1:prefixLength) ...
    - plotted_two_stage.time(1:prefixLength)));
checks.max_prefix_energy_difference = max(abs( ...
    plotted_fista.energy(1:prefixLength) ...
    - plotted_two_stage.energy(1:prefixLength)));
checks.max_prefix_pg_difference = max(abs( ...
    plotted_fista.pg(1:prefixLength) ...
    - plotted_two_stage.pg(1:prefixLength)));
checks.prefix_exactly_shared = checks.max_prefix_time_difference == 0 ...
    && checks.max_prefix_energy_difference == 0 ...
    && checks.max_prefix_pg_difference == 0;
assert(checks.prefix_exactly_shared, ...
    'The two formal curves do not share an exact common prefix.');

parameters = struct('epsilon', epsilon, 'sigma', sigma, ...
    'beta', beta, 'delta', delta, 'mass', mass, 'L', L, 'N', N, ...
    'R0', R0, 'R1', R1, ...
    'fista_plotting_cutoff', fista_plotting_cutoff);

% Standard two-stage result assembled from the exact shared FISTA prefix
% and the production Newton branch.  Figure 5.5 and Table 5.3 both consume
% this object for their common parameter case.
canonical_two_stage_result.rho = newtonFinalRho;
canonical_two_stage_result.energy = newtonFinalEnergy;
canonical_two_stage_result.target_energy = newtonFinalEnergy;
canonical_two_stage_result.problem = problem;
canonical_two_stage_result.solver = canonicalSolver;
canonical_two_stage_result.grid = problem.grid;
canonical_two_stage_result.history.main = trimHistory( ...
    mainHistory, switchIteration);
canonical_two_stage_result.history.polish = newtonHistory;
canonicalDiagnostics.main_iterations = switchIteration;
canonicalDiagnostics.main_elapsed_time = switchTime;
canonicalDiagnostics.main_energy = switchEnergy;
canonicalDiagnostics.main_pg_residual = switchPg;
canonicalDiagnostics.main_stop_reason = 'switch_to_kkt_polish';
canonicalDiagnostics.polish_iterations = newtonIterations;
canonicalDiagnostics.polish_elapsed_time = newtonElapsedTime;
canonicalDiagnostics.total_elapsed_time = switchTime + newtonElapsedTime;
canonicalDiagnostics.total_pcg_z_iterations = newtonPcgZ;
canonicalDiagnostics.total_pcg_w_iterations = newtonPcgW;
canonicalDiagnostics.max_pcg_iterations = newtonPcgMax;
canonicalDiagnostics.final_energy = newtonFinalEnergy;
canonicalDiagnostics.final_pg_residual = newtonFinalPg;
canonicalDiagnostics.final_kkt_residual = newtonFinalKkt;
canonicalDiagnostics.min_density = min(newtonFinalRho);
canonicalDiagnostics.exact_zero_count = nnz(newtonFinalRho == 0);
canonical_two_stage_result.diagnostics = canonicalDiagnostics;
canonical_timing = experiments.CanonicalSolverTiming( ...
    canonical_two_stage_result);
canonical_two_stage_result.diagnostics.canonical_solver_time = ...
    canonical_timing.total_solver_time;
canonical_two_stage_result.diagnostics.canonical_timing_protocol = ...
    canonical_timing.protocol;
canonicalResultId = ['shared_prefix_' ...
    char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'))];
canonical_two_stage_result.diagnostics.canonical_result_id = ...
    canonicalResultId;
canonical_timing.result_id = canonicalResultId;

switch_data = struct('iteration', switchIteration, ...
    'time', switchTime, 'energy', switchEnergy, 'pg', switchPg, ...
    'rho', switchState);
fista_final_pg = plotted_fista.pg(end);
two_stage_final_pg = plotted_two_stage.pg(end);
checks.figure_endpoint_time_difference = abs( ...
    plotted_two_stage.time(end) - canonical_timing.total_solver_time);
assert(checks.figure_endpoint_time_difference ...
    <= 100 * eps(max(1, canonical_timing.total_solver_time)), ...
    'The Figure 5.5 endpoint is not the canonical solver time.');

syncCanonicalEfficiency(resultFolder, repositoryRoot, ...
    canonical_two_stage_result, canonical_timing, epsilon, sigma);
table_checks = verifyCanonicalEfficiency(resultFolder, ...
    canonical_two_stage_result, canonical_timing, epsilon, sigma);
checks.table53 = table_checks;
save(dataFile, 'parameters', 'source', 'rho0', 'switch_data', ...
    'common_prefix', 'fista_continuation', 'newton_branch', ...
    'plotted_fista', 'plotted_two_stage', 'E_ref', ...
    'energyPlotFloor', 'checks', 'fista_final_pg', ...
    'two_stage_final_pg', 'canonical_two_stage_result', ...
    'canonical_timing', '-v7.3');

plotFormalHistory(plotted_fista, plotted_two_stage, switch_data, ...
    E_ref, energyPlotFloor, energyBase, residualBase);
plotIterationDiagnostic(common_prefix, diagnosticBase);
syncManuscriptFigures(energyBase, residualBase, ...
    manuscriptFigureFolder);

fprintf('\nShared-prefix solver-performance summary\n');
fprintf('  source                 : %s\n', source);
fprintf('  switch iteration       : %d\n', switchIteration);
fprintf('  switch time            : %.9f s\n', switchTime);
fprintf('  common-prefix length   : %d points (%d FISTA steps)\n', ...
    prefixLength, switchIteration);
fprintf('  two-stage final PG     : %.9e\n', two_stage_final_pg);
fprintf('  canonical solver time  : %.9f s\n', ...
    canonical_timing.total_solver_time);
fprintf('  FISTA cutoff final PG  : %.9e\n', fista_final_pg);
fprintf('  shared prefix exact    : %s\n', yesNo(checks.prefix_exactly_shared));
fprintf('  shared MAT             : %s\n', dataFile);
fprintf('  energy figure base     : %s\n', energyBase);
fprintf('  residual figure base   : %s\n', residualBase);

function config = makeConfig(beta, delta, mass, L, N, R0, R1, ...
    epsilon, sigma, r, dr, d2r, p, dp, d2p)
config = experiments.DefaultConfig();
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.mass = mass;
config.parameters.L = L;
config.parameters.N = N;
config.parameters.epsilon = epsilon;
config.fisher_regularization.r_epsilon = r;
config.fisher_regularization.dr_epsilon = dr;
config.fisher_regularization.d2r_epsilon = d2r;
config.fisher_regularization.epsilon = epsilon;
config.fisher_regularization.label = 'r_epsilon(rho)=rho+epsilon';
config.fisher_regularization.name = 'shift';
config.fisher_regularization.regularity = inf;
config.trapping_potential.V = @(x) ...
    model.HarmonicCInfPeriodicPotential(x, L, R0, R1);
config.trapping_potential.label = ...
    'harmonic core + C-infinity far-field continuation';
config.trapping_potential.boundary_periodicized = true;
config.trapping_potential.modification_start = R0;
config.trapping_potential.transition_end = R1;
config.trapping_potential.reference_V = @(x) 0.5 * x .^ 2;
config.potential_regularization.name = 'inline_p_sigma';
config.potential_regularization.label = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma';
config.potential_regularization.sigma = sigma;
config.potential_regularization.p_sigma = p;
config.potential_regularization.dp_sigma = dp;
config.potential_regularization.d2p_sigma = d2p;
config.potential_regularization.prox_type = 'generic_convex';
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.display = false;
config.solver.pg_tol = 1e-5;
config.solver.final_pg_tol = 1e-12;
config.solver.certification_tol = 1e-12;
config.solver.max_iter = 200000;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = 1e-5;
config.solver.switch.max_main_iter = 20000;
config.solver.switch.forced_pg_tol = 1e-4;
config.solver.switch.stop_at_energy_handoff = true;
config.solver.history.capture_handoff_state = true;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = 1e-12;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
end

function options = makePolishOptions(solver)
options = solver.polish;
options.projection_name = solver.projection_name;
options.projection_tol = solver.projection_tol;
options.residual_step = solver.residual_step;
options.display = solver.display;
options.entry_pg_tol = max( ...
    solver.switch.pg_entry_tol, solver.switch.forced_pg_tol);
end

function validateLegacyArchive(data, epsilon, sigma, L, N, h)
required = {'parameters', 'fistaOnlyResult', 'twoStageResult', ...
    'handoffCheck'};
assert(all(isfield(data, required)), ...
    'The legacy solver-performance archive is incomplete.');
assert(data.parameters.epsilon == epsilon ...
    && data.parameters.sigma == sigma ...
    && data.parameters.L == L && data.parameters.N == N, ...
    'The legacy archive metadata are incompatible.');
assert(data.handoffCheck.iteration_equal ...
    && data.handoffCheck.max_energy_difference == 0 ...
    && data.handoffCheck.max_pg_difference == 0, ...
    'The legacy branch prefixes were not identical.');
difference = data.fistaOnlyResult.diagnostics.handoff_rho ...
    - data.twoStageResult.diagnostics.handoff_rho;
assert(sqrt(h * sum(difference .^ 2)) == 0, ...
    'The archived Newton and FISTA branches do not share the same state.');
end

function plotFormalHistory(fista, twoStage, switchData, ERef, ...
    energyFloor, energyBase, residualBase)
plotOneHistory(fista.time, fista.energy_error_plot, ...
    twoStage.time, twoStage.energy_error_plot, switchData.time, ...
    max(abs(switchData.energy - ERef), energyFloor), ...
    '$|E-E_*|$', 'Energy convergence', energyBase);
plotOneHistory(fista.time, fista.pg, twoStage.time, twoStage.pg, ...
    switchData.time, switchData.pg, ...
    'Projected-gradient residual', 'Stationarity convergence', ...
    residualBase);
end

function plotOneHistory(fistaTime, fistaValue, twoTime, twoValue, ...
    switchTime, switchValue, yText, titleText, outputBase)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
f = figure('Color', 'w', 'Position', [100, 100, 900, 650], ...
    'Visible', 'on');
ax = axes('Parent', f);
hold(ax, 'on');
set(ax, 'YScale', 'log');
colors = colororder(ax);
semilogy(ax, fistaTime, fistaValue, '-o', ...
    'Color', colors(1, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(fistaTime)), ...
    'DisplayName', 'FISTA');
semilogy(ax, twoTime, twoValue, '--s', ...
    'Color', colors(2, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(twoTime)), ...
    'DisplayName', 'FISTA--Newton');
semilogy(ax, switchTime, switchValue, 'd', ...
    'Color', colors(3, :), 'MarkerFaceColor', colors(3, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'switch');
set(ax, 'FontSize', fontSize, 'LineWidth', 1.2);
xlabel(ax, 'Wall-clock time (s)', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
ylabel(ax, yText, 'Interpreter', 'latex', 'FontSize', fontSize);
title(ax, titleText, 'Interpreter', 'latex', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(ax, 'Location', 'best', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
xlim(ax, [0, 900]);
grid(ax, 'on');
box(ax, 'on');
saveFigureFormats(f, outputBase);
close(f);
end

function indices = markerLocations(count)
indices = unique(round(linspace(1, count, min(12, count))));
end

function saveFigureFormats(f, outputBase)
figFile = [outputBase '.fig'];
epsFile = [outputBase '.eps'];
pdfFile = [outputBase '.pdf'];
pngFile = [outputBase '.png'];
savefig(f, figFile);
set(f, 'Renderer', 'painters', 'PaperPositionMode', 'auto');
print(f, epsFile, '-depsc2', '-painters');
exportgraphics(f, pdfFile, 'ContentType', 'vector');
exportgraphics(f, pngFile, 'Resolution', 300);
end

function plotIterationDiagnostic(prefix, outputBase)
f = figure('Color', 'w', 'Position', [100, 100, 1200, 500], ...
    'Visible', 'off');
layout = tiledlayout(f, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
semilogy(ax, prefix.iteration, prefix.energy, '--', ...
    'LineWidth', 1.8, 'DisplayName', 'FISTA branch');
hold(ax, 'on');
semilogy(ax, prefix.iteration, prefix.energy, '-', ...
    'LineWidth', 1.2, 'DisplayName', 'Newton branch');
xlabel(ax, 'FISTA iteration'); ylabel(ax, 'Energy');
title(ax, 'Shared prefix: energy'); grid(ax, 'on'); box(ax, 'on');
legend(ax, 'Location', 'best');
ax = nexttile(layout);
semilogy(ax, prefix.iteration, prefix.pg, '--', ...
    'LineWidth', 1.8, 'DisplayName', 'FISTA branch');
hold(ax, 'on');
semilogy(ax, prefix.iteration, prefix.pg, '-', ...
    'LineWidth', 1.2, 'DisplayName', 'Newton branch');
xlabel(ax, 'FISTA iteration'); ylabel(ax, 'PG residual');
title(ax, 'Shared prefix: stationarity'); grid(ax, 'on'); box(ax, 'on');
legend(ax, 'Location', 'best');
savefig(f, [outputBase '.fig']);
exportgraphics(f, [outputBase '.png'], 'Resolution', 180);
close(f);
end

function syncManuscriptFigures(energyBase, residualBase, destination)
if ~isfolder(destination)
    return;
end
bases = {energyBase, residualBase};
extensions = {'.eps', '.pdf', '.png'};
for baseIndex = 1:numel(bases)
    [~, name] = fileparts(bases{baseIndex});
    for extensionIndex = 1:numel(extensions)
        extension = extensions{extensionIndex};
        copyfile([bases{baseIndex} extension], ...
            fullfile(destination, [name extension]));
    end
end
end

function history = trimHistory(history, count)
fields = fieldnames(history);
for index = 1:numel(fields)
    history.(fields{index}) = history.(fields{index})(1:count, :);
end
end

function syncCanonicalEfficiency(folder, repositoryRoot, result, timing, ...
        epsilon, sigma)
matFile = fullfile(folder, 'regularization_scale_efficiency.mat');
csvFile = fullfile(folder, 'regularization_scale_efficiency.csv');
rowFile = fullfile(folder, 'solver_regularization_cost_rows.tex');
checkpointFile = fullfile(folder, ...
    'regularization_scale_efficiency_checkpoint.mat');
if ~isfile(matFile)
    warning('SharedPrefix:MissingEfficiencyArchive', ...
        'Table 5.3 archive was not found; no table data were synchronized.');
    return;
end
data = load(matFile);
assert(all(isfield(data, {'metadata', 'records', 'rho0'})), ...
    'The Table 5.3 archive is incomplete.');
records = ensureLimitField(data.records);
index = find([records.epsilon] == epsilon & [records.sigma] == sigma);
assert(isscalar(index), 'The Table 5.3 common case is not unique.');
records(index) = updateCanonicalRecord(records(index), result, timing);
efficiencyTable = struct2table(records);
metadata = data.metadata;
metadata.timing_protocol = timing.protocol;
metadata.timing_definition = timing.definition;
metadata.common_case_source = ...
    'two_stage_vs_fista_shared_prefix.mat/canonical_two_stage_result';
rho0 = data.rho0;
% Store the very same result object consumed by Figure 5.5.  This makes
% the common Table 5.3 row traceable to one canonical data object rather
% than to a second run whose numbers merely happen to agree.
common_case_result = result;
common_case_timing = timing;
save(matFile, 'metadata', 'records', 'efficiencyTable', 'rho0', ...
    'common_case_result', 'common_case_timing', '-v7.3');
writetable(efficiencyTable, csvFile);
writeRows(rowFile, records);

manuscriptRows = fullfile(repositoryRoot, 'manuscript', 'tables', ...
    'solver_regularization_cost_rows.tex');
if isfolder(fileparts(manuscriptRows))
    copyfile(rowFile, manuscriptRows);
end

if isfile(checkpointFile)
    checkpoint = load(checkpointFile);
    if isfield(checkpoint, 'records')
        checkpoint.records = ensureLimitField(checkpoint.records);
        checkIndex = find([checkpoint.records.epsilon] == epsilon ...
            & [checkpoint.records.sigma] == sigma);
        if isscalar(checkIndex)
            checkpoint.records(checkIndex) = records(index);
            save(checkpointFile, '-struct', 'checkpoint', '-v7.3');
        end
    end
end
end

function checks = verifyCanonicalEfficiency(folder, result, timing, ...
        epsilon, sigma)
matFile = fullfile(folder, 'regularization_scale_efficiency.mat');
data = load(matFile, 'records', 'common_case_result', ...
    'common_case_timing', 'metadata');
required = {'records', 'common_case_result', 'common_case_timing', ...
    'metadata'};
assert(all(isfield(data, required)), ...
    'The canonical Table 5.3 archive is incomplete.');
index = find([data.records.epsilon] == epsilon ...
    & [data.records.sigma] == sigma);
assert(isscalar(index), 'The canonical Table 5.3 row is not unique.');
record = data.records(index);
d = result.diagnostics;

checks.same_result_id = strcmp( ...
    data.common_case_result.diagnostics.canonical_result_id, ...
    result.diagnostics.canonical_result_id);
checks.same_state = isequaln(data.common_case_result.rho, result.rho);
checks.same_main_history = isequaln( ...
    data.common_case_result.history.main, result.history.main);
checks.same_polish_history = isequaln( ...
    data.common_case_result.history.polish, result.history.polish);
checks.same_timing_struct = isequaln(data.common_case_timing, timing);
checks.time_difference = abs(record.total_solver_time ...
    - timing.total_solver_time);
checks.fista_equal = record.fista_iterations == d.main_iterations;
checks.newton_equal = record.newton_iterations == d.polish_iterations;
checks.pcg_equal = record.pcg_iterations_total == ...
    d.total_pcg_z_iterations + d.total_pcg_w_iterations;
checks.pg_difference = abs(record.final_pg - d.final_pg_residual);
checks.energy_difference = abs(record.final_energy - result.target_energy);
checks.source_equal = strcmp(data.metadata.common_case_source, ...
    'two_stage_vs_fista_shared_prefix.mat/canonical_two_stage_result');

timeTolerance = 100 * eps(max(1, timing.total_solver_time));
valueTolerance = 100 * eps(max(1, abs(result.target_energy)));
assert(checks.same_result_id && checks.same_state ...
    && checks.same_main_history && checks.same_polish_history ...
    && checks.same_timing_struct ...
    && checks.time_difference <= timeTolerance ...
    && checks.fista_equal && checks.newton_equal && checks.pcg_equal ...
    && checks.pg_difference == 0 ...
    && checks.energy_difference <= valueTolerance ...
    && checks.source_equal, ...
    ['Figure 5.5 and Table 5.3 do not consume the same canonical ' ...
    'two-stage result.']);
end

function records = ensureLimitField(records)
if ~isfield(records, 'fista_iteration_limit_reached')
    [records.fista_iteration_limit_reached] = deal(false);
end
for index = 1:numel(records)
    records(index).fista_iteration_limit_reached = ...
        records(index).fista_iterations >= 20000;
end
end

function record = updateCanonicalRecord(record, result, timing)
d = result.diagnostics;
record.fista_iterations = d.main_iterations;
record.handoff_pg = d.main_pg_residual;
record.handoff_energy = d.main_energy;
record.handoff_time = timing.fista_time;
record.newton_iterations = d.polish_iterations;
record.pcg_iterations_total = d.total_pcg_z_iterations ...
    + d.total_pcg_w_iterations;
record.pcg_iterations_mean = record.pcg_iterations_total ...
    / max(1, 2 * d.polish_iterations);
record.pcg_iterations_max = d.max_pcg_iterations;
record.fista_time = timing.fista_time;
record.newton_time = timing.newton_time;
record.total_solver_time = timing.total_solver_time;
record.final_energy = result.target_energy;
record.final_pg = d.final_pg_residual;
record.final_kkt = d.final_kkt_residual;
record.min_rho = min(result.rho);
record.exact_zero_count = nnz(result.rho == 0);
record.convergence_target_reached = ...
    d.final_pg_residual <= result.solver.final_pg_tol;
record.fista_iteration_limit_reached = ...
    d.main_iterations >= result.solver.switch.max_main_iter;
record.prox_outer_iterations_total = sum( ...
    result.history.main.prox_lambda_iterations);
record.prox_inner_iterations_total = sum( ...
    result.history.main.prox_inner_iterations);
record.timing_repeats = 1;
end

function writeRows(filename, records)
fileId = fopen(filename, 'w');
if fileId < 0
    error('Unable to open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fileId));
for index = 1:numel(records)
    record = records(index);
    fistaText = sprintf('%d', record.fista_iterations);
    if record.fista_iteration_limit_reached
        fistaText = sprintf('$%d^*$', record.fista_iterations);
    end
    rowArguments = {latexScientific(record.epsilon), ...
        latexScientific(record.sigma), fistaText, ...
        record.newton_iterations, record.pcg_iterations_total, ...
        record.total_solver_time, latexScientific(record.final_pg)};
    if index < numel(records)
        fprintf(fileId, ...
            '$%s$ & $%s$ & %s & %d & %d & %.2f & $%s$ \\\\\n', ...
            rowArguments{:});
    else
        % The manuscript supplies the final row break after \input.
        fprintf(fileId, ...
            '$%s$ & $%s$ & %s & %d & %d & %.2f & $%s$\n', ...
            rowArguments{:});
    end
    if index == 3
        fprintf(fileId, '\\midrule\n');
    end
end
end

function text = latexScientific(value)
exponent = floor(log10(abs(value)));
mantissa = value / 10 ^ exponent;
if abs(mantissa - 1) <= 10 * eps
    text = sprintf('10^{%d}', exponent);
else
    text = sprintf('%.2f\\times10^{%d}', mantissa, exponent);
end
end

function text = yesNo(value)
if value
    text = 'yes';
else
    text = 'no';
end
end
