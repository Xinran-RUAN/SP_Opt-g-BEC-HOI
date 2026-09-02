%DIAGNOSE_FIGURE55_PRE_SWITCH Diagnose the apparent pre-switch separation.
%
% This postprocessing script reads the archived production histories.  It
% does not rerun or modify either optimization trajectory.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
dataFile = fullfile(projectRoot, 'results', 'solver_performance', ...
    'two_stage_vs_fista_only.mat');
outputFolder = fullfile(projectRoot, 'results', 'solver_performance');
data = load(dataFile);

two = data.twoStageResult;
only = data.fistaOnlyResult;
twoHistory = two.history.main;
onlyHistory = only.history.main;

switchIterationTwo = two.diagnostics.handoff_iteration;
switchIterationOnly = only.diagnostics.handoff_iteration;
switchIteration = min(switchIterationTwo, switchIterationOnly);
assert(switchIterationTwo == switchIterationOnly, ...
    'The archived runs detected different switch iterations.');

twoIndex = find(twoHistory.iteration <= switchIteration);
onlyIndex = find(onlyHistory.iteration <= switchIteration);
assert(isequal(twoHistory.iteration(twoIndex), ...
    onlyHistory.iteration(onlyIndex)), ...
    'The archived FISTA iteration indices cannot be aligned.');

% Include iteration zero for the two principal trajectory comparisons.
twoCombined = data.two_stage.history.combined;
onlyCombined = data.fista_only.history.combined;
countWithInitial = switchIteration + 1;
energyTwo = twoCombined.energy(1:countWithInitial);
energyOnly = onlyCombined.energy(1:countWithInitial);
pgTwo = twoCombined.pg(1:countWithInitial);
pgOnly = onlyCombined.pg(1:countWithInitial);

energyDifference = abs(energyOnly - energyTwo);
pgDifference = abs(pgOnly - pgTwo);
energyScale = max(realmin, max(abs(energyTwo)));
pgScale = max(realmin, max(abs(pgTwo)));

metrics.switch_iteration = switchIteration;
metrics.switch_time_two_stage = two.diagnostics.handoff_time;
metrics.switch_time_fista_only = only.diagnostics.handoff_time;
metrics.figure_switch_marker_time = data.handoff_time;
metrics.switch_marker_bookkeeping_offset = abs( ...
    metrics.figure_switch_marker_time - metrics.switch_time_two_stage);
metrics.max_energy_difference = max(energyDifference);
metrics.relative_energy_difference = ...
    metrics.max_energy_difference / energyScale;
metrics.max_pg_difference = max(pgDifference);
metrics.relative_pg_difference = metrics.max_pg_difference / pgScale;

energyTolerance = 100 * eps(max(1, energyScale));
pgTolerance = 100 * eps(max(1, pgScale));
firstEnergy = find(energyDifference > energyTolerance, 1, 'first');
firstPg = find(pgDifference > pgTolerance, 1, 'first');
metrics.first_energy_divergence_iteration = indexToIteration(firstEnergy);
metrics.first_pg_divergence_iteration = indexToIteration(firstPg);

% Compare all recorded decisions that can affect the next FISTA state.
numericFields = {'accepted_L', 'tau', 'backtracks', ...
    'prox_lambda_iterations', 'prox_inner_iterations'};
for fieldIndex = 1:numel(numericFields)
    name = numericFields{fieldIndex};
    difference = abs(twoHistory.(name)(twoIndex) ...
        - onlyHistory.(name)(onlyIndex));
    metrics.(['max_' name '_difference']) = max(difference);
    first = find(difference > 0, 1, 'first');
    metrics.(['first_' name '_divergence_iteration']) = ...
        indexToIterationNoInitial(first);
end
logicalFields = {'restart', 'residual_checked', ...
    'energy_plateau_detected', 'handoff_detected'};
for fieldIndex = 1:numel(logicalFields)
    name = logicalFields{fieldIndex};
    difference = xor(twoHistory.(name)(twoIndex), ...
        onlyHistory.(name)(onlyIndex));
    metrics.([name '_mismatch_count']) = nnz(difference);
    first = find(difference, 1, 'first');
    metrics.(['first_' name '_divergence_iteration']) = ...
        indexToIterationNoInitial(first);
end

% Full state history was intentionally not archived.  The switch state was.
metrics.full_state_history_available = ...
    isfield(twoHistory, 'rho') && isfield(onlyHistory, 'rho');
if metrics.full_state_history_available
    stateDifference = twoHistory.rho(:, twoIndex) ...
        - onlyHistory.rho(:, onlyIndex);
    h = two.problem.grid.h;
    metrics.max_state_L2_difference = max(sqrt(h * ...
        sum(stateDifference .^ 2, 1)));
else
    metrics.max_state_L2_difference = NaN;
end
switchStateDifference = two.diagnostics.handoff_rho ...
    - only.diagnostics.handoff_rho;
metrics.switch_state_L2_difference = sqrt(two.problem.grid.h * ...
    sum(switchStateDifference .^ 2));

% The formal figure uses independently accumulated wall-clock histories.
timeTwo = twoHistory.elapsed_time(twoIndex);
timeOnly = onlyHistory.elapsed_time(onlyIndex);
timeDifference = abs(timeOnly - timeTwo);
relativeTimeDifference = timeDifference ./ ...
    max(max(abs(timeOnly), abs(timeTwo)), realmin);
[metrics.max_time_difference, maxTimeIndex] = max(timeDifference);
metrics.max_time_difference_iteration = ...
    twoHistory.iteration(twoIndex(maxTimeIndex));
metrics.median_relative_time_difference = median(relativeTimeDifference);
metrics.max_relative_time_difference = max(relativeTimeDifference);
metrics.switch_time_difference = abs( ...
    metrics.switch_time_fista_only - metrics.switch_time_two_stage);
metrics.switch_relative_time_difference = metrics.switch_time_difference ...
    / max(metrics.switch_time_fista_only, ...
        metrics.switch_time_two_stage);

% Verify that both plotted energy-error arrays use the single archived E_*.
expectedTwoError = abs(twoCombined.energy - data.E_star);
expectedOnlyError = abs(onlyCombined.energy - data.E_star);
metrics.two_stage_energy_reference_mismatch = max(abs( ...
    expectedTwoError - twoCombined.energy_error_raw));
metrics.fista_only_energy_reference_mismatch = max(abs( ...
    expectedOnlyError - onlyCombined.energy_error_raw));
metrics.common_energy_reference = ...
    metrics.two_stage_energy_reference_mismatch == 0 ...
    && metrics.fista_only_energy_reference_mismatch == 0;
metrics.E_star = data.E_star;

% The problems and all pre-switch algorithmic options should be identical.
metrics.problem_structs_equal_raw = isequaln(two.problem, only.problem);
metrics.problem_definitions_equal = semanticProblemEqual( ...
    two.problem, only.problem, data.rho0);
solverFields = {'name', 'splitting', 'projection_name', 'L0', 'a', ...
    'backtrack_factor', 'max_backtracks', 'feasibility_tol', ...
    'residual_check_interval'};
switchFields = {'enabled', 'energy_window', 'energy_tol', ...
    'consecutive_windows', 'min_iter', 'pg_entry_tol', ...
    'forced_pg_tol'};
metrics.common_fista_options_equal = compareFields( ...
    two.solver, only.solver, solverFields) ...
    && compareFields(two.solver.switch, only.solver.switch, switchFields);

algorithmicDifferences = [energyDifference(:); pgDifference(:); ...
    abs(twoHistory.accepted_L(twoIndex) ...
        - onlyHistory.accepted_L(onlyIndex)); ...
    abs(twoHistory.backtracks(twoIndex) ...
        - onlyHistory.backtracks(onlyIndex)); ...
    double(xor(twoHistory.restart(twoIndex), ...
        onlyHistory.restart(onlyIndex)))];
metrics.first_algorithmic_divergence_iteration = NaN;
if any(algorithmicDifferences ~= 0)
    candidateIterations = [find(energyDifference ~= 0) - 1; ...
        find(pgDifference ~= 0) - 1; ...
        find(twoHistory.accepted_L(twoIndex) ...
            ~= onlyHistory.accepted_L(onlyIndex)); ...
        find(twoHistory.backtracks(twoIndex) ...
            ~= onlyHistory.backtracks(onlyIndex)); ...
        find(twoHistory.restart(twoIndex) ...
            ~= onlyHistory.restart(onlyIndex))];
    metrics.first_algorithmic_divergence_iteration = ...
        min(candidateIterations);
end

metrics.classification = 'CASE A';
if ~metrics.common_energy_reference
    metrics.classification = 'CASE B';
elseif ~isnan(metrics.first_algorithmic_divergence_iteration) ...
        || ~metrics.common_fista_options_equal ...
        || ~metrics.problem_definitions_equal
    metrics.classification = 'CASE C';
end

% Temporary diagnostic only: the formal Fig. 5.5 remains time-based.
diagnosticFigure = figure('Color', 'w', 'Visible', 'off', ...
    'Position', [100, 100, 1200, 500]);
layout = tiledlayout(diagnosticFigure, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
iteration = (0:switchIteration).';
energyFloor = data.energyPlotFloor;
ax = nexttile(layout);
semilogy(ax, iteration, max(abs(energyOnly - data.E_star), energyFloor), ...
    '--', 'LineWidth', 1.8, 'DisplayName', 'FISTA only');
hold(ax, 'on');
semilogy(ax, iteration, max(abs(energyTwo - data.E_star), energyFloor), ...
    '-', 'LineWidth', 1.4, 'DisplayName', 'FISTA--Newton');
xlabel(ax, 'FISTA iteration'); ylabel(ax, '|E-E_*|');
title(ax, 'Pre-switch energy history'); grid(ax, 'on'); box(ax, 'on');
legend(ax, 'Location', 'best');

ax = nexttile(layout);
semilogy(ax, iteration, pgOnly, '--', 'LineWidth', 1.8, ...
    'DisplayName', 'FISTA only');
hold(ax, 'on');
semilogy(ax, iteration, pgTwo, '-', 'LineWidth', 1.4, ...
    'DisplayName', 'FISTA--Newton');
xlabel(ax, 'FISTA iteration'); ylabel(ax, 'Projected-gradient residual');
title(ax, 'Pre-switch stationarity history');
grid(ax, 'on'); box(ax, 'on'); legend(ax, 'Location', 'best');

figureFile = fullfile(outputFolder, ...
    'figure55_pre_switch_iteration_alignment.fig');
pngFile = fullfile(outputFolder, ...
    'figure55_pre_switch_iteration_alignment.png');
savefig(diagnosticFigure, figureFile);
exportgraphics(diagnosticFigure, pngFile, 'Resolution', 180);
close(diagnosticFigure);

metricsFile = fullfile(outputFolder, ...
    'figure55_pre_switch_diagnostic.mat');
save(metricsFile, 'metrics', 'energyDifference', 'pgDifference', ...
    'timeDifference', 'relativeTimeDifference');

fprintf('\nFigure 5.5 pre-switch diagnostic\n');
fprintf('  switch iteration                 : %d\n', ...
    metrics.switch_iteration);
fprintf('  switch time (two-stage)          : %.9f s\n', ...
    metrics.switch_time_two_stage);
fprintf('  switch time (FISTA only)         : %.9f s\n', ...
    metrics.switch_time_fista_only);
fprintf('  formal-figure switch marker time : %.9f s\n', ...
    metrics.figure_switch_marker_time);
fprintf('  max |Delta E|                    : %.17e\n', ...
    metrics.max_energy_difference);
fprintf('  max |Delta PG|                   : %.17e\n', ...
    metrics.max_pg_difference);
fprintf('  switch-state L2 difference       : %.17e\n', ...
    metrics.switch_state_L2_difference);
fprintf('  max accepted-L difference        : %.17e\n', ...
    metrics.max_accepted_L_difference);
fprintf('  max backtrack-count difference   : %.17e\n', ...
    metrics.max_backtracks_difference);
fprintf('  restart mismatch count           : %d\n', ...
    metrics.restart_mismatch_count);
fprintf('  max wall-clock mismatch          : %.9f s (iteration %d)\n', ...
    metrics.max_time_difference, metrics.max_time_difference_iteration);
fprintf('  switch wall-clock mismatch       : %.9f s (%.3f%%)\n', ...
    metrics.switch_time_difference, ...
    100 * metrics.switch_relative_time_difference);
fprintf('  median relative timing mismatch  : %.3f%%\n', ...
    100 * metrics.median_relative_time_difference);
fprintf('  common E_*                       : %s\n', ...
    yesNo(metrics.common_energy_reference));
fprintf('  common FISTA options             : %s\n', ...
    yesNo(metrics.common_fista_options_equal));
fprintf('  identical problem definitions    : %s\n', ...
    yesNo(metrics.problem_definitions_equal));
fprintf('  raw struct identity (handles)     : %s\n', ...
    yesNo(metrics.problem_structs_equal_raw));
fprintf('  first algorithmic divergence     : %s\n', ...
    iterationText(metrics.first_algorithmic_divergence_iteration));
fprintf('  classification                   : %s\n', ...
    metrics.classification);
fprintf('  diagnostic figure                : %s\n', pngFile);
fprintf('  diagnostic data                  : %s\n', metricsFile);

function value = indexToIteration(index)
if isempty(index)
    value = NaN;
else
    value = index - 1;
end
end

function value = indexToIterationNoInitial(index)
if isempty(index)
    value = NaN;
else
    value = index;
end
end

function tf = compareFields(left, right, names)
tf = true;
for fieldIndex = 1:numel(names)
    name = names{fieldIndex};
    tf = tf && isfield(left, name) && isfield(right, name) ...
        && isequaln(left.(name), right.(name));
end
end

function tf = semanticProblemEqual(left, right, rho0)
tf = isequaln(left.grid, right.grid) ...
    && isequaln(left.V, right.V) ...
    && isequaln(left.beta, right.beta) ...
    && isequaln(left.delta, right.delta) ...
    && isequaln(left.mass, right.mass) ...
    && isequaln(left.entropy, right.entropy) ...
    && isequaln(left.fisher_regularization.epsilon, ...
        right.fisher_regularization.epsilon) ...
    && isequaln(left.potential_regularization.sigma, ...
        right.potential_regularization.sigma);
testRho = [0; rho0(:); ...
    left.fisher_regularization.epsilon; ...
    left.potential_regularization.sigma];
fisherNames = {'s_epsilon', 'ds_epsilon', 'd2s_epsilon'};
potentialNames = {'p_sigma', 'dp_sigma', 'd2p_sigma'};
for fieldIndex = 1:numel(fisherNames)
    name = fisherNames{fieldIndex};
    tf = tf && isequaln( ...
        left.fisher_regularization.(name)(testRho), ...
        right.fisher_regularization.(name)(testRho));
end
for fieldIndex = 1:numel(potentialNames)
    name = potentialNames{fieldIndex};
    tf = tf && isequaln( ...
        left.potential_regularization.(name)(testRho), ...
        right.potential_regularization.(name)(testRho));
end
end

function text = yesNo(tf)
if tf
    text = 'yes';
else
    text = 'no';
end
end

function text = iterationText(value)
if isnan(value)
    text = 'none';
else
    text = sprintf('%d', value);
end
end
