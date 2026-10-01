%RUN_SOLVER_SPLITTING_COMPARISON Add smooth-potential FISTA to Figure 5.5.
%
% All three trajectories minimize the same E_{epsilon,sigma,N}.  The
% production curves are loaded from the shared-prefix archive.  The only
% algorithmic difference in the new branch is whether V*p_sigma belongs
% to the smooth term or to the proximal term.

if ~exist('reuse_smooth_result', 'var')
    reuse_smooth_result = true;
end
clearvars -except reuse_smooth_result;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
startup_HOI();

sharedFile = fullfile(projectRoot, 'results', 'solver_performance', ...
    'two_stage_vs_fista_shared_prefix.mat');
dataFile = fullfile(projectRoot, 'results', ...
    'solver_performance_with_smooth_potential.mat');
figureFolder = fullfile(projectRoot, 'figs');
repositoryRoot = fileparts(projectRoot);
manuscriptFigureFolder = fullfile(repositoryRoot, 'manuscript', 'figs');
energyBase = fullfile(figureFolder, ...
    'solver_energy_history_shared_prefix');
residualBase = fullfile(figureFolder, ...
    'solver_residual_history_shared_prefix');
maxIterations = 50000;
plotTimeWindow = [0, 900];

assert(isfile(sharedFile), ...
    'Run RUN_TWO_STAGE_SOLVER_PERFORMANCE_SHARED_PREFIX first.');
shared = load(sharedFile);
validateSharedData(shared);
parameters = shared.parameters;
assert(parameters.epsilon == 1e-2 && parameters.sigma == 1e-4 ...
    && parameters.L == 32 && parameters.N == 4096, ...
    'The shared-prefix archive does not have the Figure 5.5 parameters.');

problem = shared.canonical_two_stage_result.problem;
rho0 = shared.rho0;
smoothSolver = shared.canonical_two_stage_result.solver;
smoothSolver.name = 'fista_cd';
smoothSolver.splitting = 'smooth_potential';
smoothSolver.polish_mode = 'none';
smoothSolver.switch.enabled = true;
smoothSolver.switch.stop_at_energy_handoff = false;
smoothSolver.max_iter = maxIterations;
smoothSolver.display = false;

% Same-objective validation.  The smooth splitting reconstructs the
% production target from SmoothEnergy plus the potential contribution;
% its analytic gradient must equal the production full gradient.
testRho = 0.65 * rho0 + 0.35 * shared.newton_branch.final_rho;
testRho = testRho * (problem.mass ...
    / src.constraints.Mass(testRho, problem.grid.h));
fullEnergy = src.discretization.ps.Energy(testRho, problem);
kineticSmoothEnergy = src.discretization.ps.SmoothEnergy(testRho, problem);
[potentialDensity, potentialDerivative] = src.potential.Evaluate( ...
    testRho, problem.potential_regularization, ...
    problem.fisher_regularization.epsilon);
reconstructedEnergy = kineticSmoothEnergy + problem.grid.h * sum( ...
    problem.V(:) .* potentialDensity(:));
fullGradient = src.discretization.ps.Gradient(testRho, problem);
reconstructedGradient = src.discretization.ps.SmoothGradient( ...
    testRho, problem) + problem.V(:) .* potentialDerivative(:);
validation.energy_objective_difference = abs( ...
    fullEnergy - reconstructedEnergy);
validation.gradient_infinity_difference = max(abs( ...
    fullGradient - reconstructedGradient));
validation.gradient_relative_L2_difference = norm( ...
    fullGradient - reconstructedGradient) / max(1, norm(fullGradient));
assert(validation.energy_objective_difference <= ...
    500 * eps(max(1, abs(fullEnergy))), ...
    'The two splittings do not reconstruct the same objective.');
assert(validation.gradient_infinity_difference <= ...
    500 * eps(max(1, norm(fullGradient, inf))), ...
    'The smooth-potential gradient is inconsistent with the target.');

loadedSmooth = false;
if reuse_smooth_result && isfile(dataFile)
    old = load(dataFile, 'smooth_result', 'parameters', 'maxIterations');
    if isfield(old, 'smooth_result') && isfield(old, 'parameters') ...
            && isfield(old, 'maxIterations') ...
            && compatibleParameters(old.parameters, parameters) ...
            && old.maxIterations == maxIterations
        smoothResult = old.smooth_result;
        loadedSmooth = true;
        fprintf('Reused compatible smooth-potential FISTA result.\n');
    end
end
if ~loadedSmooth
    fprintf(['Running smooth-potential FISTA: epsilon=%.1e, sigma=%.1e, ' ...
        'L=%g, N=%d, maxit=%d ...\n'], parameters.epsilon, ...
        parameters.sigma, parameters.L, parameters.N, maxIterations);
    smoothResult = src.SolveGroundState1D(problem, rho0, smoothSolver);
    % Save the expensive branch immediately; postprocessing below can be
    % repeated without rerunning the solver.
    save(dataFile, 'parameters', 'maxIterations', 'smoothResult', ...
        'validation', '-v7.3');
end

smoothHistory = smoothResult.history.main;
initialEnergy = src.discretization.ps.Energy(rho0, problem);
initialGradient = src.discretization.ps.Gradient(rho0, problem);
initialPg = src.solvers.FullGradientMapping( ...
    rho0, initialGradient, problem, smoothSolver);
smooth.time = [0; smoothHistory.elapsed_time(:)];
smooth.energy = [initialEnergy; smoothHistory.augmented_energy(:)];
smooth.pg = [initialPg; smoothHistory.full_pg_residual(:)];
smooth.iteration = [0; smoothHistory.iteration(:)];

E_ref = shared.E_ref;
energyPlotFloor = max(1e-16, 10 * eps(max(1, abs(E_ref))));
smooth.energy_error_raw = abs(smooth.energy - E_ref);
smooth.energy_error_plot = max(smooth.energy_error_raw, energyPlotFloor);

proxPrefixL = shared.common_prefix.accepted_L(2:end);
proxPrefixBacktracks = shared.common_prefix.backtracks(2:end);
smoothL = smoothHistory.accepted_L(:);
smoothBacktracks = smoothHistory.backtracks(:);
step_statistics.prox_prefix.median_tau = median(1 ./ proxPrefixL);
step_statistics.prox_prefix.min_tau = min(1 ./ proxPrefixL);
step_statistics.prox_prefix.max_tau = max(1 ./ proxPrefixL);
step_statistics.prox_prefix.total_backtracking = sum(proxPrefixBacktracks);
step_statistics.smooth.median_tau = median(1 ./ smoothL);
step_statistics.smooth.min_tau = min(1 ./ smoothL);
step_statistics.smooth.max_tau = max(1 ./ smoothL);
step_statistics.smooth.total_backtracking = sum(smoothBacktracks);

final_comparison.smooth_energy_error = abs( ...
    smoothResult.target_energy - E_ref);
final_comparison.smooth_state_L2_difference = sqrt(problem.grid.h * sum( ...
    (smoothResult.rho - shared.newton_branch.final_rho) .^ 2));
final_comparison.smooth_final_pg = ...
    smoothResult.diagnostics.final_pg_residual;
final_comparison.smooth_final_iteration = ...
    smoothResult.diagnostics.main_iterations;
final_comparison.smooth_solver_time = ...
    smoothResult.diagnostics.main_elapsed_time;
final_comparison.smooth_mass_error = ...
    smoothResult.diagnostics.mass_error;
final_comparison.smooth_min_rho = min(smoothResult.rho);

plotThreeHistories(shared.plotted_fista, shared.plotted_two_stage, ...
    smooth, shared.switch_data, E_ref, energyPlotFloor, plotTimeWindow, ...
    energyBase, residualBase);
if isfolder(manuscriptFigureFolder)
    syncRenderedFigures({energyBase, residualBase}, ...
        manuscriptFigureFolder);
end

source_shared_prefix_file = sharedFile; %#ok<NASGU>
prox_fista = shared.plotted_fista; %#ok<NASGU>
prox_two_stage = shared.plotted_two_stage; %#ok<NASGU>
switch_data = shared.switch_data; %#ok<NASGU>
save(dataFile, 'parameters', 'maxIterations', 'plotTimeWindow', ...
    'source_shared_prefix_file', 'rho0', 'E_ref', 'energyPlotFloor', ...
    'prox_fista', 'prox_two_stage', 'switch_data', 'smooth', ...
    'smoothResult', 'step_statistics', 'validation', ...
    'final_comparison', '-v7.3');

fprintf('\nSolver-splitting comparison summary\n');
fprintf('  same-objective energy defect : %.3e\n', ...
    validation.energy_objective_difference);
fprintf('  full-gradient infinity defect: %.3e\n', ...
    validation.gradient_infinity_difference);
fprintf('  smooth FISTA iterations      : %d\n', ...
    final_comparison.smooth_final_iteration);
fprintf('  smooth FISTA solver time     : %.9f s\n', ...
    final_comparison.smooth_solver_time);
fprintf('  smooth FISTA final full PG   : %.9e\n', ...
    final_comparison.smooth_final_pg);
fprintf('  smooth FISTA final dE        : %.9e\n', ...
    final_comparison.smooth_energy_error);
fprintf('  smooth/reference state L2    : %.9e\n', ...
    final_comparison.smooth_state_L2_difference);
fprintf('  prox-prefix tau median/min   : %.3e / %.3e\n', ...
    step_statistics.prox_prefix.median_tau, ...
    step_statistics.prox_prefix.min_tau);
fprintf('  prox-prefix backtracks       : %d\n', ...
    step_statistics.prox_prefix.total_backtracking);
fprintf('  smooth tau median/min        : %.3e / %.3e\n', ...
    step_statistics.smooth.median_tau, step_statistics.smooth.min_tau);
fprintf('  smooth backtracks            : %d\n', ...
    step_statistics.smooth.total_backtracking);
fprintf('  data MAT                     : %s\n', dataFile);
fprintf('  energy figure base           : %s\n', energyBase);
fprintf('  residual figure base         : %s\n', residualBase);

function validateSharedData(data)
required = {'parameters', 'rho0', 'E_ref', 'common_prefix', ...
    'plotted_fista', 'plotted_two_stage', 'switch_data', ...
    'newton_branch', 'canonical_two_stage_result'};
assert(all(isfield(data, required)), ...
    'The shared-prefix archive is incomplete.');
assert(data.checks.prefix_exactly_shared, ...
    'The archived production trajectories do not share an exact prefix.');
end

function tf = compatibleParameters(a, b)
names = {'epsilon', 'sigma', 'beta', 'delta', 'mass', 'L', 'N'};
tf = all(cellfun(@(name) isfield(a, name) && isfield(b, name) ...
    && isequal(a.(name), b.(name)), names));
end

function plotThreeHistories(proxFista, proxTwoStage, smooth, ...
    switchData, ERef, energyFloor, timeWindow, energyBase, residualBase)
plotOne(proxFista.time, proxFista.energy_error_plot, ...
    proxTwoStage.time, proxTwoStage.energy_error_plot, ...
    smooth.time, smooth.energy_error_plot, switchData.time, ...
    max(abs(switchData.energy - ERef), energyFloor), ...
    '$|E-E_*|$', 'Energy convergence', timeWindow, energyBase);
plotOne(proxFista.time, proxFista.pg, proxTwoStage.time, ...
    proxTwoStage.pg, smooth.time, smooth.pg, switchData.time, ...
    switchData.pg, 'Projected-gradient residual', ...
    'Stationarity convergence', timeWindow, residualBase);
end

function plotOne(proxTime, proxValue, twoTime, twoValue, smoothTime, ...
    smoothValue, switchTime, switchValue, yText, titleText, ...
    timeWindow, outputBase)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
f = figure('Color', 'w', 'Position', [100, 100, 900, 650], ...
    'Visible', 'on');
ax = axes('Parent', f);
hold(ax, 'on');
set(ax, 'YScale', 'log');
colors = colororder(ax);
semilogy(ax, proxTime, proxValue, '-o', 'Color', colors(1, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(proxTime)), ...
    'DisplayName', 'potential-prox FISTA');
semilogy(ax, twoTime, twoValue, '--s', 'Color', colors(2, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(twoTime)), ...
    'DisplayName', 'potential-prox FISTA--Newton');
semilogy(ax, smoothTime, smoothValue, '-.^', 'Color', colors(3, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(smoothTime)), ...
    'DisplayName', 'smooth-potential FISTA');
semilogy(ax, switchTime, switchValue, 'd', 'Color', colors(4, :), ...
    'MarkerFaceColor', colors(4, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, 'DisplayName', 'switch');
set(ax, 'FontSize', fontSize, 'LineWidth', 1.2);
xlabel(ax, 'Wall-clock time (s)', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
ylabel(ax, yText, 'Interpreter', 'latex', 'FontSize', fontSize);
title(ax, titleText, 'Interpreter', 'latex', 'FontSize', fontSize, ...
    'FontWeight', 'normal');
legend(ax, 'Location', 'best', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
xlim(ax, timeWindow);
grid(ax, 'on');
box(ax, 'on');
saveFigureFormats(f, outputBase);
close(f);
end

function indices = markerLocations(count)
indices = unique(round(linspace(1, count, min(12, count))));
end

function saveFigureFormats(f, outputBase)
savefig(f, [outputBase '.fig']);
set(f, 'Renderer', 'painters', 'PaperPositionMode', 'auto');
print(f, [outputBase '.eps'], '-depsc2', '-painters');
exportgraphics(f, [outputBase '.pdf'], 'ContentType', 'vector');
exportgraphics(f, [outputBase '.png'], 'Resolution', 300);
end

function syncRenderedFigures(bases, destination)
extensions = {'.eps', '.pdf', '.png'};
for baseIndex = 1:numel(bases)
    [~, name] = fileparts(bases{baseIndex});
    for extensionIndex = 1:numel(extensions)
        extension = extensions{extensionIndex};
        copyfile([bases{baseIndex}, extension], ...
            fullfile(destination, [name, extension]));
    end
end
end
