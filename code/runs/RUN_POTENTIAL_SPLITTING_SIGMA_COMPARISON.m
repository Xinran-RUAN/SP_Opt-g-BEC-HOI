%RUN_POTENTIAL_SPLITTING_SIGMA_COMPARISON Fixed-budget splitting study.
%
% Four first-order runs minimize the same sigma-dependent target energy.
% They stop only when the full projected-gradient residual reaches 1e-6
% or the canonical solver wall-clock reaches 900 seconds.  max_iter is a
% nonbinding API safeguard and is checked explicitly after every run.

if ~exist('reuse_splitting_results', 'var')
    reuse_splitting_results = true;
end
clearvars -except reuse_splitting_results;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
repositoryRoot = fileparts(projectRoot);
addpath(projectRoot);
startup_HOI();

epsilon = 1e-2;
sigmaList = [1e-4, 1e-6];
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 4096;
pgTolerance = 1e-6;
timeLimit = 900;
nonbindingMaxIterations = 500000;
thresholds = [1e-3, 1e-4, 1e-5, 1e-6];
splittingList = {'potential_prox', 'smooth_potential'};

resultFolder = fullfile(projectRoot, 'results');
solverResultFolder = fullfile(resultFolder, 'solver_performance');
figureFolder = fullfile(projectRoot, 'figs');
manuscriptFigureFolder = fullfile(repositoryRoot, 'manuscript', 'figs');
if ~isfolder(resultFolder), mkdir(resultFolder); end
if ~isfolder(figureFolder), mkdir(figureFolder); end

dataFile = fullfile(resultFolder, ...
    'potential_splitting_sigma_comparison.mat');
checkpointFile = fullfile(resultFolder, ...
    'potential_splitting_sigma_comparison_checkpoint.mat');
summaryCsv = fullfile(resultFolder, ...
    'potential_splitting_sigma_comparison_summary.csv');
targetCsv = fullfile(resultFolder, ...
    'potential_splitting_sigma_time_to_target.csv');
sharedFile = fullfile(solverResultFolder, ...
    'two_stage_vs_fista_shared_prefix.mat');
assert(isfile(sharedFile), ...
    'The certified Figure 5.5 shared-prefix archive is required.');
shared = load(sharedFile);
validateSharedArchive(shared, epsilon, 1e-4, beta, delta, mass, L, N);

rho0 = shared.rho0;
baseProblem = shared.canonical_two_stage_result.problem;
baseSolver = shared.canonical_two_stage_result.solver;
signature = struct('epsilon', epsilon, 'sigma_list', sigmaList, ...
    'beta', beta, 'delta', delta, 'mass', mass, 'L', L, 'N', N, ...
    'pg_tolerance', pgTolerance, 'time_limit', timeLimit, ...
    'nonbinding_max_iterations', nonbindingMaxIterations, ...
    'splitting_list', {splittingList});
checkpoint.signature = signature;
checkpoint.runs = cell(numel(sigmaList), numel(splittingList));
if reuse_splitting_results && isfile(checkpointFile)
    loaded = load(checkpointFile, 'checkpoint');
    if isfield(loaded, 'checkpoint') ...
            && isfield(loaded.checkpoint, 'signature') ...
            && isequaln(loaded.checkpoint.signature, signature)
        checkpoint = loaded.checkpoint;
        fprintf('Reused compatible fixed-budget splitting checkpoint.\n');
    else
        fprintf('Ignored incompatible splitting checkpoint.\n');
    end
end

objectiveChecks = repmat(struct( ...
    'sigma', NaN, 'energy_reconstruction_difference', NaN, ...
    'full_gradient_infinity_difference', NaN, ...
    'full_gradient_relative_L2_difference', NaN), numel(sigmaList), 1);
records = repmat(emptyRecord(numel(thresholds)), ...
    numel(sigmaList), numel(splittingList));
pairComparison = repmat(struct('sigma', NaN, ...
    'energy_difference', NaN, 'state_L2_difference', NaN), ...
    numel(sigmaList), 1);

for sigmaIndex = 1:numel(sigmaList)
    sigma = sigmaList(sigmaIndex);
    problem = problemAtSigma(baseProblem, sigma);
    representativeRho = representativeState( ...
        rho0, referenceState(sigma, shared, projectRoot, N));
    objectiveChecks(sigmaIndex) = verifyObjective( ...
        representativeRho, problem, sigma);
    assert(objectiveChecks(sigmaIndex).energy_reconstruction_difference ...
        <= 500 * eps(max(1, abs( ...
        src.discretization.ps.Energy(representativeRho, problem)))), ...
        'Energy reconstruction failed for sigma %.3e.', sigma);
    assert(objectiveChecks(sigmaIndex).full_gradient_infinity_difference ...
        <= 500 * eps(max(1, norm( ...
        src.discretization.ps.Gradient(representativeRho, problem), inf))), ...
        'Gradient reconstruction failed for sigma %.3e.', sigma);

    % Warm up target evaluation and both proximal backends outside the
    % canonical timers, avoiding first-call FFT/JIT asymmetry.
    warmupProblem(problem, rho0, baseSolver);

    for splittingIndex = 1:numel(splittingList)
        splitting = splittingList{splittingIndex};
        if isempty(checkpoint.runs{sigmaIndex, splittingIndex})
            solver = fixedBudgetSolver(baseSolver, splitting, ...
                pgTolerance, timeLimit, nonbindingMaxIterations);
            fprintf(['Running sigma=%.0e, %-18s, PG<=%.1e or ' ...
                'time>=%.0fs ...\n'], sigma, splitting, ...
                pgTolerance, timeLimit);
            result = src.SolveGroundState1D(problem, rho0, solver);
            checkpoint.runs{sigmaIndex, splittingIndex} = result;
            save(checkpointFile, 'checkpoint', '-v7.3');
        else
            result = checkpoint.runs{sigmaIndex, splittingIndex};
            fprintf('Reused sigma=%.0e, %s.\n', sigma, splitting);
        end
        records(sigmaIndex, splittingIndex) = makeRecord( ...
            result, rho0, sigma, splitting, thresholds, pgTolerance, ...
            timeLimit, nonbindingMaxIterations);
    end
end

references = cell(numel(sigmaList), 1);
for sigmaIndex = 1:numel(sigmaList)
    sigma = sigmaList(sigmaIndex);
    references{sigmaIndex} = referenceState( ...
        sigma, shared, projectRoot, N);
    for splittingIndex = 1:numel(splittingList)
        result = checkpoint.runs{sigmaIndex, splittingIndex};
        records(sigmaIndex, splittingIndex).reference_energy_error = abs( ...
            result.target_energy - references{sigmaIndex}.target_energy);
        records(sigmaIndex, splittingIndex).reference_state_L2_error = ...
            sqrt(baseProblem.grid.h * sum( ...
            (result.rho - references{sigmaIndex}.rho) .^ 2));
    end
    if all(strcmp({records(sigmaIndex, :).termination_reason}, ...
            'stationarity_reached'))
        proxResult = checkpoint.runs{sigmaIndex, 1};
        smoothResult = checkpoint.runs{sigmaIndex, 2};
        pairComparison(sigmaIndex).sigma = sigma;
        pairComparison(sigmaIndex).energy_difference = abs( ...
            proxResult.target_energy - smoothResult.target_energy);
        pairComparison(sigmaIndex).state_L2_difference = sqrt( ...
            baseProblem.grid.h * sum( ...
            (proxResult.rho - smoothResult.rho) .^ 2));
    else
        pairComparison(sigmaIndex).sigma = sigma;
        pairComparison(sigmaIndex).energy_difference = NaN;
        pairComparison(sigmaIndex).state_L2_difference = NaN;
    end
end

medianTauRatio = [records(:, 1).median_tau] ./ ...
    [records(:, 2).median_tau];
timeRatios = computeTimeRatios(records, thresholds);

plotSplittingDiagnostics(records, figureFolder, timeLimit);
plotFigure55(shared, records(1, 2), figureFolder, timeLimit);
syncFigure55(figureFolder, manuscriptFigureFolder);

[summaryTable, timeTargetTable] = makeTables(records, thresholds);
writetable(summaryTable, summaryCsv);
writetable(timeTargetTable, targetCsv);
save(dataFile, 'signature', 'rho0', 'objectiveChecks', 'records', ...
    'references', 'pairComparison', 'medianTauRatio', 'timeRatios', ...
    'summaryTable', 'timeTargetTable', 'thresholds', ...
    'checkpointFile', 'sharedFile', '-v7.3');

printSummary(records, objectiveChecks, thresholds, medianTauRatio, ...
    timeRatios, dataFile, summaryCsv, targetCsv, figureFolder, ...
    nonbindingMaxIterations);

function record = emptyRecord(thresholdCount)
record = struct('sigma', NaN, 'method', '', 'splitting', '', ...
    'termination_reason', '', 'final_time', NaN, ...
    'final_iteration', NaN, 'final_pg', NaN, 'final_energy', NaN, ...
    'final_rho', [], 'mass_error', NaN, 'min_rho', NaN, ...
    'median_tau', NaN, 'min_tau', NaN, 'max_tau', NaN, ...
    'mean_tau', NaN, 'total_backtracks', NaN, ...
    'thresholds', nan(1, thresholdCount), ...
    'threshold_iterations', nan(1, thresholdCount), ...
    'threshold_times', nan(1, thresholdCount), ...
    'time', [], 'iteration', [], 'energy', [], 'pg', [], ...
    'tau', [], 'backtracks', [], 'reference_energy_error', NaN, ...
    'reference_state_L2_error', NaN, 'max_iter_triggered', false);
end

function validateSharedArchive(data, epsilon, sigma, beta, delta, mass, L, N)
required = {'parameters', 'rho0', 'E_ref', 'plotted_fista', ...
    'plotted_two_stage', 'switch_data', 'canonical_two_stage_result'};
assert(all(isfield(data, required)), ...
    'The Figure 5.5 shared-prefix archive is incomplete.');
p = data.parameters;
assert(p.epsilon == epsilon && p.sigma == sigma && p.beta == beta ...
    && p.delta == delta && p.mass == mass && p.L == L && p.N == N, ...
    'The Figure 5.5 shared-prefix metadata are incompatible.');
end

function problem = problemAtSigma(baseProblem, sigma)
problem = baseProblem;
problem.potential_regularization.name = 'inline_p_sigma';
problem.potential_regularization.label = sprintf( ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma, sigma=%.3e', sigma);
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(rho) ...
    rho ./ hypot(rho, sigma);
problem.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
problem.potential_regularization.prox_type = 'generic_convex';
end

function rho = representativeState(rho0, reference)
rho = 0.65 * rho0 + 0.35 * reference.rho;
rho = rho * (reference.problem.mass ...
    / src.constraints.Mass(rho, reference.problem.grid.h));
end

function check = verifyObjective(rho, problem, sigma)
fullEnergy = src.discretization.ps.Energy(rho, problem);
smoothEnergy = src.discretization.ps.SmoothEnergy(rho, problem);
[potentialDensity, potentialDerivative] = src.potential.Evaluate( ...
    rho, problem.potential_regularization, ...
    problem.fisher_regularization.epsilon);
reconstructedEnergy = smoothEnergy + problem.grid.h * sum( ...
    problem.V(:) .* potentialDensity(:));
fullGradient = src.discretization.ps.Gradient(rho, problem);
reconstructedGradient = src.discretization.ps.SmoothGradient( ...
    rho, problem) + problem.V(:) .* potentialDerivative(:);
check.sigma = sigma;
check.energy_reconstruction_difference = abs( ...
    fullEnergy - reconstructedEnergy);
check.full_gradient_infinity_difference = max(abs( ...
    fullGradient - reconstructedGradient));
check.full_gradient_relative_L2_difference = norm( ...
    fullGradient - reconstructedGradient) / max(1, norm(fullGradient));
end

function warmupProblem(problem, rho0, baseSolver)
energy = src.discretization.ps.Energy(rho0, problem);
gradient = src.discretization.ps.Gradient(rho0, problem);
src.solvers.FullGradientMapping(rho0, gradient, problem, baseSolver);
src.discretization.ps.SmoothEnergy(rho0, problem);
smoothGradient = src.discretization.ps.SmoothGradient(rho0, problem);
for splitting = {'potential_prox', 'smooth_potential'}
    solver = baseSolver;
    solver.splitting = splitting{1};
    src.solvers.ProxStep(rho0, smoothGradient, 1, problem, solver);
end
assert(isfinite(energy));
end

function solver = fixedBudgetSolver(base, splitting, pgTol, timeLimit, maxIter)
solver = base;
solver.name = 'fista_cd';
solver.splitting = splitting;
solver.polish_mode = 'none';
solver.pg_tol = pgTol;
solver.final_pg_tol = pgTol;
solver.certification_tol = pgTol;
solver.max_iter = maxIter;
solver.time_limit = timeLimit;
solver.switch.enabled = true;
solver.switch.stop_at_energy_handoff = false;
solver.history.capture_handoff_state = false;
solver.display = false;
end

function record = makeRecord(result, rho0, sigma, splitting, thresholds, ...
    pgTol, timeLimit, maxIter)
history = result.history.main;
initialEnergy = src.discretization.ps.Energy(rho0, result.problem);
initialGradient = src.discretization.ps.Gradient(rho0, result.problem);
initialPg = src.solvers.FullGradientMapping( ...
    rho0, initialGradient, result.problem, result.solver);
record = emptyRecord(numel(thresholds));
record.sigma = sigma;
record.splitting = splitting;
if strcmp(splitting, 'potential_prox')
    record.method = 'prox-FISTA';
else
    record.method = 'smooth-FISTA';
end
record.final_time = result.diagnostics.main_elapsed_time;
record.final_iteration = result.diagnostics.main_iterations;
record.final_pg = result.diagnostics.final_pg_residual;
record.final_energy = result.target_energy;
record.final_rho = result.rho;
record.mass_error = result.diagnostics.mass_error;
record.min_rho = min(result.rho);
record.time = [0; history.elapsed_time(:)];
record.iteration = [0; history.iteration(:)];
record.energy = [initialEnergy; history.augmented_energy(:)];
record.pg = [initialPg; history.full_pg_residual(:)];
record.tau = history.tau(:);
record.backtracks = history.backtracks(:);
record.median_tau = median(record.tau);
record.min_tau = min(record.tau);
record.max_tau = max(record.tau);
record.mean_tau = mean(record.tau);
record.total_backtracks = sum(record.backtracks);
record.thresholds = thresholds;
for index = 1:numel(thresholds)
    first = find(record.pg <= thresholds(index), 1);
    if ~isempty(first)
        record.threshold_iterations(index) = record.iteration(first);
        record.threshold_times(index) = record.time(first);
    end
end
stop = result.diagnostics.main_stop_reason;
if record.final_pg <= pgTol && strcmp(stop, 'final_pg_tolerance')
    record.termination_reason = 'stationarity_reached';
elseif strcmp(stop, 'time_limit') && record.final_time >= timeLimit
    record.termination_reason = 'time_limit';
else
    error(['Unexpected termination for sigma %.3e, %s: %s, ' ...
        'iter=%d, PG=%.3e, time=%.3f.'], sigma, splitting, stop, ...
        record.final_iteration, record.final_pg, record.final_time);
end
record.max_iter_triggered = record.final_iteration >= maxIter;
assert(~record.max_iter_triggered, ...
    'The nonbinding max_iter safeguard unexpectedly triggered.');
end

function reference = referenceState(sigma, shared, projectRoot, N)
if sigma == 1e-4
    reference = shared.canonical_two_stage_result;
    return;
end
filename = fullfile(projectRoot, 'results', ...
    'sigma_spectral_accuracy', ...
    'sigma_mesh_1em06_eps1em02_L32_checkpoint.mat');
assert(isfile(filename), ...
    'The existing sigma=1e-6 high-accuracy mesh checkpoint is missing.');
loaded = load(filename, 'checkpoint');
solutions = loaded.checkpoint.solutions;
reference = [];
for index = 1:numel(solutions)
    if ~isempty(solutions{index}) && solutions{index}.grid.N == N
        reference = solutions{index};
        break;
    end
end
assert(~isempty(reference), ...
    'No N=%d sigma=1e-6 reference was found.', N);
assert(reference.parameters.epsilon == 1e-2 ...
    && reference.grid.L == 32 ...
    && reference.problem.potential_regularization.sigma == sigma, ...
    'The sigma=1e-6 reference metadata are incompatible.');
end

function ratios = computeTimeRatios(records, thresholds)
ratios = nan(numel(records(:, 1)), numel(thresholds));
for sigmaIndex = 1:size(records, 1)
    for thresholdIndex = 1:numel(thresholds)
        proxTime = records(sigmaIndex, 1).threshold_times(thresholdIndex);
        smoothTime = records(sigmaIndex, 2).threshold_times(thresholdIndex);
        if isfinite(proxTime) && isfinite(smoothTime)
            ratios(sigmaIndex, thresholdIndex) = smoothTime / proxTime;
        end
    end
end
end

function plotSplittingDiagnostics(records, folder, timeLimit)
styles = {'-o', '--s', '-.^', ':d'};
labels = {'$\sigma=10^{-4}$, prox-FISTA', ...
    '$\sigma=10^{-4}$, smooth-FISTA', ...
    '$\sigma=10^{-6}$, prox-FISTA', ...
    '$\sigma=10^{-6}$, smooth-FISTA'};
ordered = [records(1, 1), records(1, 2), records(2, 1), records(2, 2)];
f = newFigure();
ax = axes('Parent', f, 'YScale', 'log'); hold(ax, 'on');
colors = colororder(ax);
for index = 1:numel(ordered)
    indices = sampleIndices(numel(ordered(index).time), 3000);
    semilogy(ax, ordered(index).time(indices), ...
        ordered(index).pg(indices), styles{index}, ...
        'Color', colors(index, :), 'LineWidth', 2, 'MarkerSize', 8, ...
        'MarkerIndices', markerLocations(numel(indices)), ...
        'DisplayName', labels{index});
end
formatAxis(ax, 'Wall-clock time (s)', ...
    'Projected-gradient residual', ...
    'Potential-splitting stationarity comparison');
xlim(ax, [0, timeLimit]);
legend(ax, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
saveFigure(f, fullfile(folder, 'splitting_sigma_comparison_pg'));
close(f);

f = newFigure();
ax = axes('Parent', f, 'YScale', 'log'); hold(ax, 'on');
colors = colororder(ax);
for index = 1:numel(ordered)
    tauTime = ordered(index).time(2:end);
    indices = sampleIndices(numel(tauTime), 3000);
    semilogy(ax, tauTime(indices), ordered(index).tau(indices), ...
        styles{index}, 'Color', colors(index, :), 'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'MarkerIndices', markerLocations(numel(indices)), ...
        'DisplayName', labels{index});
end
formatAxis(ax, 'Wall-clock time (s)', 'Accepted step size $\tau_k$', ...
    'Accepted FISTA step sizes');
xlim(ax, [0, timeLimit]);
legend(ax, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
saveFigure(f, fullfile(folder, 'splitting_sigma_comparison_tau'));
close(f);
end

function plotFigure55(shared, smoothRecord, folder, timeLimit)
Eref = shared.E_ref;
floorValue = max(1e-16, 10 * eps(max(1, abs(Eref))));
smoothEnergyError = max(abs(smoothRecord.energy - Eref), floorValue);
plotFigure55One(shared.plotted_fista.time, ...
    shared.plotted_fista.energy_error_plot, ...
    shared.plotted_two_stage.time, ...
    shared.plotted_two_stage.energy_error_plot, ...
    smoothRecord.time, smoothEnergyError, shared.switch_data.time, ...
    max(abs(shared.switch_data.energy - Eref), floorValue), ...
    '$|E-E_*|$', 'Energy convergence', timeLimit, ...
    fullfile(folder, 'solver_energy_history_shared_prefix'));
plotFigure55One(shared.plotted_fista.time, shared.plotted_fista.pg, ...
    shared.plotted_two_stage.time, shared.plotted_two_stage.pg, ...
    smoothRecord.time, smoothRecord.pg, shared.switch_data.time, ...
    shared.switch_data.pg, 'Projected-gradient residual', ...
    'Stationarity convergence', timeLimit, ...
    fullfile(folder, 'solver_residual_history_shared_prefix'));
end

function plotFigure55One(proxTime, proxValue, twoTime, twoValue, ...
    smoothTime, smoothValue, switchTime, switchValue, yText, ...
    titleText, timeLimit, outputBase)
f = newFigure();
ax = axes('Parent', f, 'YScale', 'log'); hold(ax, 'on');
colors = colororder(ax);
semilogy(ax, proxTime, proxValue, '-o', 'Color', colors(1, :), ...
    'LineWidth', 2, 'MarkerSize', 9, ...
    'MarkerIndices', markerLocations(numel(proxTime)), ...
    'DisplayName', 'potential-prox FISTA');
semilogy(ax, twoTime, twoValue, '--s', 'Color', colors(2, :), ...
    'LineWidth', 2, 'MarkerSize', 9, ...
    'MarkerIndices', markerLocations(numel(twoTime)), ...
    'DisplayName', 'potential-prox FISTA--Newton');
smoothIndices = sampleIndices(numel(smoothTime), 5000);
semilogy(ax, smoothTime(smoothIndices), smoothValue(smoothIndices), ...
    '-.^', 'Color', colors(3, :), 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerIndices', markerLocations(numel(smoothIndices)), ...
    'DisplayName', 'smooth-potential FISTA');
semilogy(ax, switchTime, switchValue, 'd', 'Color', colors(4, :), ...
    'MarkerFaceColor', colors(4, :), 'LineWidth', 2, ...
    'MarkerSize', 9, 'DisplayName', 'switch');
formatAxis(ax, 'Wall-clock time (s)', yText, titleText);
xlim(ax, [0, timeLimit]);
legend(ax, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
saveFigure(f, outputBase);
close(f);
end

function f = newFigure()
f = figure('Color', 'w', 'Position', [100, 100, 900, 650], ...
    'Visible', 'on');
end

function formatAxis(ax, xText, yText, titleText)
set(ax, 'FontSize', 20, 'LineWidth', 1.2);
xlabel(ax, xText, 'Interpreter', 'latex', 'FontSize', 20);
ylabel(ax, yText, 'Interpreter', 'latex', 'FontSize', 20);
title(ax, titleText, 'Interpreter', 'latex', 'FontSize', 20, ...
    'FontWeight', 'normal');
grid(ax, 'on'); box(ax, 'on');
end

function saveFigure(f, base)
savefig(f, [base, '.fig']);
set(f, 'Renderer', 'painters', 'PaperPositionMode', 'auto');
print(f, [base, '.eps'], '-depsc2', '-painters');
exportgraphics(f, [base, '.png'], 'Resolution', 300);
end

function indices = sampleIndices(count, maximum)
indices = unique(round(linspace(1, count, min(count, maximum))));
end

function indices = markerLocations(count)
indices = unique(round(linspace(1, count, min(12, count))));
end

function syncFigure55(source, destination)
if ~isfolder(destination), return; end
names = {'solver_energy_history_shared_prefix', ...
    'solver_residual_history_shared_prefix'};
extensions = {'.eps', '.png'};
for nameIndex = 1:numel(names)
    for extensionIndex = 1:numel(extensions)
        extension = extensions{extensionIndex};
        copyfile(fullfile(source, [names{nameIndex}, extension]), ...
            fullfile(destination, [names{nameIndex}, extension]));
    end
end
end

function [summary, targets] = makeTables(records, thresholds)
flat = [records(1, 1), records(1, 2), records(2, 1), records(2, 2)];
flat = flat(:);
summary = table([flat.sigma].', string({flat.method}).', ...
    string({flat.termination_reason}).', [flat.final_time].', ...
    [flat.final_iteration].', [flat.final_pg].', [flat.median_tau].', ...
    [flat.min_tau].', [flat.max_tau].', [flat.mean_tau].', ...
    [flat.total_backtracks].', 'VariableNames', ...
    {'sigma', 'method', 'termination', 'final_time', 'final_iteration', ...
    'final_pg', 'median_tau', 'min_tau', 'max_tau', 'mean_tau', ...
    'total_backtracks'});
targetMatrix = vertcat(flat.threshold_times);
targets = table([flat.sigma].', string({flat.method}).', ...
    targetMatrix(:, 1), targetMatrix(:, 2), targetMatrix(:, 3), ...
    targetMatrix(:, 4), 'VariableNames', {'sigma', 'method', ...
    'T_pg_1em3', 'T_pg_1em4', 'T_pg_1em5', 'T_pg_1em6'});
assert(isequal(thresholds, [1e-3, 1e-4, 1e-5, 1e-6]));
end

function printSummary(records, checks, thresholds, tauRatios, timeRatios, ...
    dataFile, summaryCsv, targetCsv, figureFolder, maxIter)
fprintf('\nFixed-budget potential-splitting comparison\n');
fprintf(['sigma method       termination           time(s)  iter       PG' ...
    '       median tau  min tau    max tau    mean tau   backtracks\n']);
for sigmaIndex = 1:size(records, 1)
    for methodIndex = 1:size(records, 2)
        r = records(sigmaIndex, methodIndex);
        fprintf(['%.0e %-12s %-20s %8.3f %7d %.3e %.3e %.3e ' ...
            '%.3e %.3e %d\n'], r.sigma, r.method, ...
            r.termination_reason, r.final_time, r.final_iteration, ...
            r.final_pg, r.median_tau, r.min_tau, r.max_tau, ...
            r.mean_tau, r.total_backtracks);
    end
end
fprintf('\nTime to full-PG thresholds (s; -- means not reached)\n');
fprintf('sigma method       ');
fprintf(' PG<=%.0e', thresholds);
fprintf('\n');
for sigmaIndex = 1:size(records, 1)
    for methodIndex = 1:size(records, 2)
        r = records(sigmaIndex, methodIndex);
        fprintf('%.0e %-12s', r.sigma, r.method);
        for index = 1:numel(thresholds)
            if isfinite(r.threshold_times(index))
                fprintf(' %10.3f', r.threshold_times(index));
            else
                fprintf(' %10s', '--');
            end
        end
        fprintf('\n');
    end
end
for sigmaIndex = 1:size(records, 1)
    fprintf(['sigma=%.0e: p_sigma''''(0)=%.1e, median tau ratio ' ...
        'prox/smooth=%.3e\n'], records(sigmaIndex, 1).sigma, ...
        1 / records(sigmaIndex, 1).sigma, tauRatios(sigmaIndex));
    fprintf('  objective defects: energy %.3e, gradient inf %.3e\n', ...
        checks(sigmaIndex).energy_reconstruction_difference, ...
        checks(sigmaIndex).full_gradient_infinity_difference);
    fprintf('  time ratios smooth/prox:');
    fprintf(' %.3g', timeRatios(sigmaIndex, :));
    fprintf('\n');
end
fprintf('Nonbinding max_iter safeguard: %d (not triggered).\n', maxIter);
fprintf('Data MAT  : %s\n', dataFile);
fprintf('Summary CSV: %s\n', summaryCsv);
fprintf('Target CSV : %s\n', targetCsv);
fprintf('Figure dir : %s\n', figureFolder);
end
