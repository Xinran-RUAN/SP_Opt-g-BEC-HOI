%RUN_TWO_STAGE_SOLVER_PERFORMANCE Two-stage solver performance study.
%
% Experiment A compares production FISTA->Newton-PCG with the identical
% FISTA iteration continued beyond the energy-window switch. Experiment B
% measures the two-stage cost at five fixed-N regularization-scale cases.
% No continuation warm starts are used: every timed solve starts from rho0.

if ~exist('performance_smoke_test', 'var')
    performance_smoke_test = false;
end
if ~exist('timing_repeats', 'var')
    timing_repeats = 1;
end
clearvars -except performance_smoke_test timing_repeats;
clc;

% =========================== editable setup ===========================
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 4096;

epsilon_A = 1e-2;
sigma_A = 1e-4;
epsilon_list = [1e-2; 1e-3; 1e-4];
sigma_list = [1e-2; 1e-3; 1e-4];

R0_fraction = 0.75;
R1_fraction = 0.90;
fista_only_max_iter = 50000;
reuse_completed_results = true;

% Simple Fisher formula remains directly editable in this run file.
r_epsilon = @(rho, epsilon) rho + epsilon;
dr_epsilon = @(rho, epsilon) ones(size(rho));
d2r_epsilon = @(rho, epsilon) zeros(size(rho));

% Stable potential smoothing formulas.
p_sigma = @(rho, sigma) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho, sigma) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho, sigma) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
% =====================================================================

if performance_smoke_test
    N = 512;
    fista_only_max_iter = 5000;
end
if ~isscalar(timing_repeats) || timing_repeats < 1 ...
        || timing_repeats ~= round(timing_repeats)
    error('timing_repeats must be a positive integer.');
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
resultFolder = fullfile(root, 'results', 'solver_performance');
if performance_smoke_test
    resultFolder = fullfile(resultFolder, 'smoke');
    figureFolder = fullfile(resultFolder, 'figs');
    fileSuffix = '_smoke';
else
    figureFolder = fullfile(root, 'figs');
    fileSuffix = '';
end
if ~isfolder(resultFolder), mkdir(resultFolder); end
if ~isfolder(figureFolder), mkdir(figureFolder); end

twoStageFile = fullfile(resultFolder, ...
    ['two_stage_vs_fista_only', fileSuffix, '.mat']);
efficiencyMatFile = fullfile(resultFolder, ...
    ['regularization_scale_efficiency', fileSuffix, '.mat']);
efficiencyCsvFile = fullfile(resultFolder, ...
    ['regularization_scale_efficiency', fileSuffix, '.csv']);
latexRowsFile = fullfile(resultFolder, ...
    ['solver_regularization_cost_rows', fileSuffix, '.tex']);
efficiencyCheckpointFile = fullfile(resultFolder, ...
    ['regularization_scale_efficiency_checkpoint', fileSuffix, '.mat']);
sharedPrefixFile = fullfile(resultFolder, ...
    ['two_stage_vs_fista_shared_prefix', fileSuffix, '.mat']);
energyFigFile = fullfile(figureFolder, ...
    ['solver_energy_history', fileSuffix, '.fig']);
energyEpsFile = fullfile(figureFolder, ...
    ['solver_energy_history', fileSuffix, '.eps']);
residualFigFile = fullfile(figureFolder, ...
    ['solver_residual_history', fileSuffix, '.fig']);
residualEpsFile = fullfile(figureFolder, ...
    ['solver_residual_history', fileSuffix, '.eps']);

R0 = R0_fraction * L;
R1 = R1_fraction * L;
formula.r_epsilon = r_epsilon;
formula.dr_epsilon = dr_epsilon;
formula.d2r_epsilon = d2r_epsilon;
formula.p_sigma = p_sigma;
formula.dp_sigma = dp_sigma;
formula.d2p_sigma = d2p_sigma;

baseConfig = makeConfig(beta, delta, mass, L, N, R0, R1, ...
    epsilon_A, sigma_A, formula);
grid = model.SetupGrid1D(baseConfig.parameters);
rho0 = model.InitialDensity(grid, mass, 'semismooth');
assert(min(rho0) >= 0, 'The common initial density is not nonnegative.');
assert(abs(src.constraints.Mass(rho0, grid.h) - mass) <= 1e-12, ...
    'The common initial density does not have the requested mass.');

% Untimed warm-up: compile the actual energy/gradient/FFT/PG/prox/PCG path.
fprintf('Untimed production warm-up ...\n');
warmConfig = makeConfig(beta, delta, mass, L, 64, ...
    R0, R1, epsilon_A, sigma_A, formula);
warmGrid = model.SetupGrid1D(warmConfig.parameters);
warmRho = model.InitialDensity(warmGrid, mass, 'semismooth');
warmConfig.solver.display = false;
experiments.SolveGroundState(warmConfig, warmRho);

%% Experiment A: FISTA-only versus production two-stage.
fprintf('\nEXPERIMENT A: two-stage solve (epsilon %.1e, sigma %.1e, N %d)\n', ...
    epsilon_A, sigma_A, N);
twoStageConfig = baseConfig;
twoStageConfig.solver.history.capture_handoff_state = true;
resumeA = false;
if reuse_completed_results && isfile(twoStageFile)
    savedA = load(twoStageFile, 'parameters', 'twoStageResult', ...
        'twoStageTiming', 'fistaOnlyResult', 'fistaOnlyTiming');
    resumeA = all(isfield(savedA, {'parameters', 'twoStageResult', ...
        'twoStageTiming', 'fistaOnlyResult', 'fistaOnlyTiming'})) ...
        && savedA.parameters.N == N ...
        && savedA.parameters.L == L ...
        && savedA.parameters.epsilon == epsilon_A ...
        && savedA.parameters.sigma == sigma_A;
end
if resumeA
    fprintf('Reusing compatible completed Experiment A.\n');
    twoStageResult = savedA.twoStageResult;
    twoStageTiming = savedA.twoStageTiming;
else
    [twoStageResult, twoStageTiming] = timedSolve( ...
        twoStageConfig, rho0, timing_repeats);
end

productionPgTarget = twoStageConfig.solver.final_pg_tol;
pgTwoStageFinal = twoStageResult.diagnostics.final_pg_residual;
pgTargetCompare = max(productionPgTarget, 1.05 * pgTwoStageFinal);

fistaOnlyConfig = twoStageConfig;
fistaOnlyConfig.solver.polish_mode = 'none';
fistaOnlyConfig.solver.switch.stop_at_energy_handoff = false;
fistaOnlyConfig.solver.max_iter = fista_only_max_iter;
fistaOnlyConfig.solver.final_pg_tol = pgTargetCompare;
fistaOnlyConfig.solver.certification_tol = pgTargetCompare;
fprintf(['EXPERIMENT A: FISTA-only diagnostic, continuing past switch ' ...
    'to PG %.3e or %d iterations\n'], ...
    pgTargetCompare, fista_only_max_iter);
if resumeA
    fistaOnlyResult = savedA.fistaOnlyResult;
    fistaOnlyTiming = savedA.fistaOnlyTiming;
else
    [fistaOnlyResult, fistaOnlyTiming] = timedSolve( ...
        fistaOnlyConfig, rho0, timing_repeats);
end

handoffCheck = compareHandoff(twoStageResult, fistaOnlyResult, grid.h);
assert(handoffCheck.iteration_equal, ...
    'FISTA-only and two-stage detected different switch iterations.');
assert(handoffCheck.max_energy_difference <= 100 * eps( ...
    max(1, abs(twoStageResult.diagnostics.main_energy))), ...
    'Pre-switch energy histories differ by %.3e.', ...
    handoffCheck.max_energy_difference);
assert(handoffCheck.max_pg_difference <= 1e-13, ...
    'Pre-switch PG histories differ by %.3e.', ...
    handoffCheck.max_pg_difference);
assert(handoffCheck.state_L2_difference <= 1e-13, ...
    'Switch states differ by %.3e in L2.', ...
    handoffCheck.state_L2_difference);

E_star = twoStageResult.target_energy;
twoStageTrajectory = buildTrajectory(twoStageResult, rho0);
fistaOnlyTrajectory = buildTrajectory(fistaOnlyResult, rho0);
energyPlotFloor = max(1e-16, 10 * eps(max(1, abs(E_star))));
twoStageTrajectory.energy_error_raw = abs( ...
    twoStageTrajectory.energy - E_star);
twoStageTrajectory.energy_error_plot = max( ...
    twoStageTrajectory.energy_error_raw, energyPlotFloor);
fistaOnlyTrajectory.energy_error_raw = abs( ...
    fistaOnlyTrajectory.energy - E_star);
fistaOnlyTrajectory.energy_error_plot = max( ...
    fistaOnlyTrajectory.energy_error_raw, energyPlotFloor);

handoff_time = twoStageResult.diagnostics.main_elapsed_time;
handoff_pg = twoStageResult.diagnostics.main_pg_residual;
handoff_energy = twoStageResult.diagnostics.main_energy;
handoff_fista_iter = twoStageResult.diagnostics.main_iterations;
fista_only_reached_common_pg = ...
    fistaOnlyResult.diagnostics.final_pg_residual <= pgTargetCompare;
two_stage_total_time = twoStageResult.diagnostics.total_elapsed_time;
fista_only_total_time = fistaOnlyResult.diagnostics.total_elapsed_time;

two_stage.history.main = twoStageResult.history.main;
two_stage.history.polish = twoStageResult.history.polish;
two_stage.history.combined = twoStageTrajectory;
two_stage.final_state = twoStageResult.rho;
two_stage.diagnostics = twoStageResult.diagnostics;
two_stage.timing = twoStageTiming;
fista_only.history.main = fistaOnlyResult.history.main;
fista_only.history.polish = fistaOnlyResult.history.polish;
fista_only.history.combined = fistaOnlyTrajectory;
fista_only.final_state = fistaOnlyResult.rho;
fista_only.diagnostics = fistaOnlyResult.diagnostics;
fista_only.timing = fistaOnlyTiming;

parameters = struct('epsilon', epsilon_A, 'sigma', sigma_A, ...
    'beta', beta, 'delta', delta, 'mass', mass, 'L', L, 'N', N, ...
    'R0', R0, 'R1', R1, 'timing_repeats', timing_repeats, ...
    'smoke_test', performance_smoke_test, ...
    'production_pg_target', productionPgTarget, ...
    'pg_target_compare', pgTargetCompare);
pg_two_stage_final = pgTwoStageFinal;
pg_fista_only_final = fistaOnlyResult.diagnostics.final_pg_residual;
save(twoStageFile, 'parameters', 'rho0', 'E_star', 'two_stage', ...
    'fista_only', 'handoff_time', 'handoff_pg', 'handoff_energy', ...
    'handoff_fista_iter', 'pg_two_stage_final', ...
    'pg_fista_only_final', 'fista_only_reached_common_pg', ...
    'two_stage_total_time', 'fista_only_total_time', ...
    'handoffCheck', 'energyPlotFloor', 'twoStageResult', ...
    'twoStageTiming', 'fistaOnlyResult', 'fistaOnlyTiming', '-v7.3');

plotHistories(twoStageTrajectory, fistaOnlyTrajectory, ...
    handoff_time, handoff_energy, handoff_pg, E_star, ...
    energyPlotFloor, energyFigFile, energyEpsFile, ...
    residualFigFile, residualEpsFile);
reportSanity('two-stage', twoStageResult, mass);
reportSanity('FISTA-only', fistaOnlyResult, mass);
reportEnergyIncrease('two-stage FISTA', ...
    twoStageResult.history.main.augmented_energy);
reportEnergyIncrease('FISTA-only', ...
    fistaOnlyResult.history.main.augmented_energy);
if isfield(twoStageResult.history.polish, 'energy')
    reportEnergyIncrease('Newton polish', ...
        twoStageResult.history.polish.energy);
end

fprintf('\nRepresentative comparison\n');
fprintf('  switch FISTA iter  : %d\n', handoff_fista_iter);
fprintf('  switch PG          : %.6e\n', handoff_pg);
fprintf('  Newton iter        : %d\n', ...
    twoStageResult.diagnostics.polish_iterations);
fprintf('  total PCG          : %d\n', totalPcg(twoStageResult));
fprintf('  two-stage final PG : %.6e\n', pgTwoStageFinal);
fprintf('  two-stage time     : %.6f s\n', two_stage_total_time);
fprintf('  FISTA-only iter    : %d\n', ...
    fistaOnlyResult.diagnostics.main_iterations);
fprintf('  FISTA-only final PG: %.6e\n', pg_fista_only_final);
fprintf('  FISTA-only time    : %.6f s\n', fista_only_total_time);
fprintf('  common PG reached  : %s\n', yesNo(fista_only_reached_common_pg));
fprintf('  switch state L2    : %.3e\n', ...
    handoffCheck.state_L2_difference);

%% Experiment B: fixed-N regularization-scale efficiency.
% N is deliberately fixed: this isolates solver difficulty from the cost
% of changing the number of Fourier degrees of freedom.
caseEpsilon = [1e-2; 1e-3; 1e-4; 1e-2; 1e-2];
caseSigma = [1e-2; 1e-2; 1e-2; 1e-3; 1e-4];
records = repmat(emptyRecord(), numel(caseEpsilon), 1);
completedCases = false(numel(caseEpsilon), 1);
efficiencySignature = struct('N', N, 'L', L, 'beta', beta, ...
    'delta', delta, 'mass', mass, 'R0', R0, 'R1', R1, ...
    'case_epsilon', caseEpsilon, 'case_sigma', caseSigma, ...
    'timing_repeats', timing_repeats);
if reuse_completed_results && isfile(efficiencyCheckpointFile)
    checkpoint = load(efficiencyCheckpointFile, 'signature', ...
        'records', 'completedCases');
    if all(isfield(checkpoint, {'signature', 'records', ...
            'completedCases'})) ...
            && isequaln(checkpoint.signature, efficiencySignature) ...
            && numel(checkpoint.records) == numel(records)
        records = checkpoint.records;
        completedCases = checkpoint.completedCases;
        fprintf('Reusing %d completed efficiency cases.\n', ...
            nnz(completedCases));
    end
end
if ~isfield(records, 'fista_iteration_limit_reached')
    [records.fista_iteration_limit_reached] = deal(false);
end
for caseIndex = 1:numel(caseEpsilon)
    epsilon = caseEpsilon(caseIndex);
    sigma = caseSigma(caseIndex);
    canonicalCommon = epsilon == epsilon_A && sigma == sigma_A ...
        && isfile(sharedPrefixFile);
    if canonicalCommon
        shared = load(sharedPrefixFile, 'parameters', ...
            'canonical_two_stage_result', 'canonical_timing');
        assert(isfield(shared, 'canonical_two_stage_result') ...
            && isfield(shared, 'canonical_timing'), ...
            'The shared-prefix archive predates canonical timing data.');
        assert(shared.parameters.epsilon == epsilon ...
            && shared.parameters.sigma == sigma ...
            && shared.parameters.L == L && shared.parameters.N == N, ...
            'The shared-prefix common case has incompatible metadata.');
        state = shared.canonical_two_stage_result;
        timing = shared.canonical_timing;
        fprintf(['\nUsing the shared-prefix canonical result for ' ...
            'epsilon %.1e, sigma %.1e.\n'], epsilon, sigma);
    elseif completedCases(caseIndex)
        fprintf('\nReused EXPERIMENT B case %d/%d: epsilon %.1e, sigma %.1e\n', ...
            caseIndex, numel(caseEpsilon), epsilon, sigma);
        continue;
    elseif epsilon == epsilon_A && sigma == sigma_A
        state = twoStageResult;
        timing = twoStageTiming;
        fprintf('\nReusing Experiment A for epsilon %.1e, sigma %.1e.\n', ...
            epsilon, sigma);
    else
        config = makeConfig(beta, delta, mass, L, N, R0, R1, ...
            epsilon, sigma, formula);
        fprintf('\nEXPERIMENT B case %d/%d: epsilon %.1e, sigma %.1e\n', ...
            caseIndex, numel(caseEpsilon), epsilon, sigma);
        [state, timing] = timedSolve(config, rho0, timing_repeats);
    end
    records(caseIndex) = summarizeCase(state, timing, epsilon, sigma);
    completedCases(caseIndex) = true;
    signature = efficiencySignature;
    save(efficiencyCheckpointFile, 'signature', 'records', ...
        'completedCases', '-v7.3');
    reportSanity(sprintf('case %d', caseIndex), state, mass);
end

efficiencyTable = struct2table(records);
writetable(efficiencyTable, efficiencyCsvFile);
metadata = parameters;
metadata.description = ['Fixed-N two-stage efficiency; every case ' ...
    'starts from the same projected Gaussian rho0.'];
metadata.timing_protocol = 'canonical_solver_time_v1';
metadata.timing_definition = [ ...
    'production FISTA plus production KKT refinement; setup, I/O, ' ...
    'plotting, and experiment-only postprocessing excluded'];
metadata.fisher = 'r_epsilon(rho)=rho+epsilon';
metadata.potential = ...
    'p_sigma(rho)=rho^2/(hypot(rho,sigma)+sigma)';
save(efficiencyMatFile, 'metadata', 'records', ...
    'efficiencyTable', 'rho0', '-v7.3');
writeLatexRows(latexRowsFile, records);
printEfficiencyTable(records);

fprintf('\nSaved solver-performance outputs\n');
fprintf('  %s\n', twoStageFile);
fprintf('  %s\n', efficiencyMatFile);
fprintf('  %s\n', efficiencyCsvFile);
fprintf('  %s\n', latexRowsFile);
fprintf('  %s\n', energyFigFile);
fprintf('  %s\n', energyEpsFile);
fprintf('  %s\n', residualFigFile);
fprintf('  %s\n', residualEpsFile);

function config = makeConfig(beta, delta, mass, L, N, R0, R1, ...
    epsilon, sigma, formula)
config = experiments.DefaultConfig();
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.mass = mass;
config.parameters.L = L;
config.parameters.N = N;
config.parameters.epsilon = epsilon;
config.fisher_regularization.r_epsilon = ...
    @(rho) formula.r_epsilon(rho, epsilon);
config.fisher_regularization.dr_epsilon = ...
    @(rho) formula.dr_epsilon(rho, epsilon);
config.fisher_regularization.d2r_epsilon = ...
    @(rho) formula.d2r_epsilon(rho, epsilon);
config.fisher_regularization.epsilon = epsilon;
config.fisher_regularization.label = ...
    'r_epsilon(rho)=rho+epsilon';
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
config.potential_regularization.p_sigma = ...
    @(rho) formula.p_sigma(rho, sigma);
config.potential_regularization.dp_sigma = ...
    @(rho) formula.dp_sigma(rho, sigma);
config.potential_regularization.d2p_sigma = ...
    @(rho) formula.d2p_sigma(rho, sigma);
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
config.solver.history.capture_handoff_state = false;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = 1e-12;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
end

function [result, timing] = timedSolve(config, rho0, repeats)
[result, timing] = experiments.CanonicalTimedSolve(config, rho0, repeats);
end

function diagnostic = compareHandoff(twoStage, fistaOnly, h)
twoDiagnostic = twoStage.diagnostics;
fistaDiagnostic = fistaOnly.diagnostics;
diagnostic.iteration_equal = twoDiagnostic.handoff_iteration ...
    == fistaDiagnostic.handoff_iteration;
count = min(twoDiagnostic.handoff_iteration, ...
    fistaDiagnostic.handoff_iteration);
twoEnergy = twoStage.history.main.augmented_energy(1:count);
fistaEnergy = fistaOnly.history.main.augmented_energy(1:count);
twoPg = twoStage.history.main.full_pg_residual(1:count);
fistaPg = fistaOnly.history.main.full_pg_residual(1:count);
diagnostic.max_energy_difference = max(abs(twoEnergy - fistaEnergy));
diagnostic.max_pg_difference = max(abs(twoPg - fistaPg));
if isempty(twoDiagnostic.handoff_rho) ...
        || isempty(fistaDiagnostic.handoff_rho)
    diagnostic.state_L2_difference = Inf;
else
    difference = twoDiagnostic.handoff_rho ...
        - fistaDiagnostic.handoff_rho;
    diagnostic.state_L2_difference = sqrt(h * sum(difference .^ 2));
end
end

function trajectory = buildTrajectory(result, rho0)
problem = result.problem;
solver = result.solver;
initialEnergy = src.discretization.ps.Energy(rho0, problem);
initialGradient = src.discretization.ps.Gradient(rho0, problem);
initialPg = src.solvers.FullGradientMapping( ...
    rho0, initialGradient, problem, solver);
main = result.history.main;
trajectory.time = [0; main.elapsed_time(:)];
trajectory.energy = [initialEnergy; main.augmented_energy(:)];
trajectory.pg = [initialPg; main.full_pg_residual(:)];
trajectory.stage = [repmat({'FISTA'}, numel(trajectory.time), 1)];
if result.diagnostics.polish_entered ...
        && isfield(result.history.polish, 'elapsed_time')
    polish = result.history.polish;
    trajectory.time = [trajectory.time; ...
        result.diagnostics.main_elapsed_time + polish.elapsed_time(:)];
    trajectory.energy = [trajectory.energy; polish.energy(:)];
    trajectory.pg = [trajectory.pg; polish.full_pg_residual(:)];
    trajectory.stage = [trajectory.stage; ...
        repmat({'Newton'}, numel(polish.elapsed_time), 1)];
end
end

function plotHistories(twoStage, fistaOnly, handoffTime, ...
    handoffEnergy, handoffPg, EStar, energyFloor, ...
    energyFig, energyEps, residualFig, residualEps)
fontSize = 20;
lineWidth = 2;
markerSize = 9;

f = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
ax = axes('Parent', f); hold(ax, 'on');
set(ax, 'YScale', 'log');
colors = colororder(ax);
semilogy(ax, fistaOnly.time, fistaOnly.energy_error_plot, '--o', ...
    'Color', colors(1, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(fistaOnly.time)), ...
    'DisplayName', 'FISTA only');
semilogy(ax, twoStage.time, twoStage.energy_error_plot, '-s', ...
    'Color', colors(2, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(twoStage.time)), ...
    'DisplayName', 'FISTA--Newton');
semilogy(ax, handoffTime, max(abs(handoffEnergy - EStar), energyFloor), ...
    'd', 'Color', colors(3, :), 'MarkerFaceColor', colors(3, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'switch');
formatHistoryAxes(ax, 'Wall-clock time (s)', '$|E-E_*|$', ...
    'Energy convergence', fontSize);
savePaperFigure(f, energyFig, energyEps);

f = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
ax = axes('Parent', f); hold(ax, 'on');
set(ax, 'YScale', 'log');
colors = colororder(ax);
semilogy(ax, fistaOnly.time, fistaOnly.pg, '--o', ...
    'Color', colors(1, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(fistaOnly.time)), ...
    'DisplayName', 'FISTA only');
semilogy(ax, twoStage.time, twoStage.pg, '-s', ...
    'Color', colors(2, :), 'LineWidth', lineWidth, ...
    'MarkerSize', markerSize, ...
    'MarkerIndices', markerLocations(numel(twoStage.time)), ...
    'DisplayName', 'FISTA--Newton');
semilogy(ax, handoffTime, handoffPg, 'd', ...
    'Color', colors(3, :), 'MarkerFaceColor', colors(3, :), ...
    'LineWidth', lineWidth, 'MarkerSize', markerSize, ...
    'DisplayName', 'switch');
formatHistoryAxes(ax, 'Wall-clock time (s)', ...
    'Projected-gradient residual', 'Stationarity convergence', fontSize);
savePaperFigure(f, residualFig, residualEps);
end

function indices = markerLocations(count)
indices = unique(round(linspace(1, count, min(12, count))));
end

function formatHistoryAxes(ax, xText, yText, titleText, fontSize)
set(ax, 'FontSize', fontSize, 'LineWidth', 1.2);
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, xText, 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(ax, yText, 'Interpreter', 'latex', 'FontSize', fontSize);
title(ax, titleText, 'Interpreter', 'latex', ...
    'FontSize', fontSize, 'FontWeight', 'normal');
legend(ax, 'Location', 'best', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
end

function savePaperFigure(f, figFile, epsFile)
savefig(f, figFile);
screenDpi = get(groot, 'ScreenPixelsPerInch');
figurePosition = get(f, 'Position');
paperSize = figurePosition(3:4) / screenDpi;
set(f, 'Renderer', 'painters', 'PaperUnits', 'inches', ...
    'PaperSize', paperSize, 'PaperPosition', [0, 0, paperSize], ...
    'PaperPositionMode', 'manual');
print(f, epsFile, '-depsc2', '-painters');
close(f);
end

function record = emptyRecord()
record = struct('epsilon', NaN, 'sigma', NaN, ...
    'fista_iterations', NaN, 'handoff_pg', NaN, ...
    'handoff_energy', NaN, 'handoff_time', NaN, ...
    'newton_iterations', NaN, 'pcg_iterations_total', NaN, ...
    'pcg_iterations_mean', NaN, 'pcg_iterations_max', NaN, ...
    'fista_time', NaN, 'newton_time', NaN, ...
    'total_solver_time', NaN, 'final_energy', NaN, ...
    'final_pg', NaN, 'final_kkt', NaN, 'min_rho', NaN, ...
    'exact_zero_count', NaN, 'convergence_target_reached', false, ...
    'fista_iteration_limit_reached', false, ...
    'prox_outer_iterations_total', NaN, ...
    'prox_inner_iterations_total', NaN, ...
    'timing_repeats', NaN);
end

function record = summarizeCase(result, timing, epsilon, sigma)
d = result.diagnostics;
record = emptyRecord();
record.epsilon = epsilon;
record.sigma = sigma;
record.fista_iterations = d.main_iterations;
record.handoff_pg = d.main_pg_residual;
record.handoff_energy = d.main_energy;
record.handoff_time = d.main_elapsed_time;
record.newton_iterations = d.polish_iterations;
record.pcg_iterations_total = totalPcg(result);
if d.polish_iterations > 0
    record.pcg_iterations_mean = record.pcg_iterations_total ...
        / (2 * d.polish_iterations);
else
    record.pcg_iterations_mean = 0;
end
record.pcg_iterations_max = d.max_pcg_iterations;
record.fista_time = d.main_elapsed_time;
record.newton_time = d.polish_elapsed_time;
record.total_solver_time = d.total_elapsed_time;
record.final_energy = result.target_energy;
record.final_pg = d.final_pg_residual;
record.final_kkt = d.final_kkt_residual;
record.min_rho = min(result.rho);
record.exact_zero_count = nnz(result.rho == 0);
record.convergence_target_reached = d.final_pg_residual ...
    <= result.solver.final_pg_tol;
record.fista_iteration_limit_reached = d.main_iterations ...
    >= result.solver.switch.max_main_iter;
record.prox_outer_iterations_total = sum( ...
    result.history.main.prox_lambda_iterations);
record.prox_inner_iterations_total = sum( ...
    result.history.main.prox_inner_iterations);
record.timing_repeats = timing.repeats;
end

function value = totalPcg(result)
value = result.diagnostics.total_pcg_z_iterations ...
    + result.diagnostics.total_pcg_w_iterations;
end

function printEfficiencyTable(records)
fprintf('\n-----------------------------------------------------------------------\n');
fprintf(' epsilon     sigma      FISTA  Newton   PCG     time(s)      final PG\n');
fprintf('-----------------------------------------------------------------------\n');
for j = 1:numel(records)
    r = records(j);
    fprintf(' %.1e   %.1e   %6d  %6d  %6d  %10.4f   %.3e\n', ...
        r.epsilon, r.sigma, r.fista_iterations, ...
        r.newton_iterations, r.pcg_iterations_total, ...
        r.total_solver_time, r.final_pg);
end
fprintf('-----------------------------------------------------------------------\n');
end

function writeLatexRows(file, records)
fileId = fopen(file, 'w');
if fileId < 0
    error('Unable to open %s for writing.', file);
end
cleanup = onCleanup(@() fclose(fileId));
for j = 1:numel(records)
    r = records(j);
    fistaText = sprintf('%d', r.fista_iterations);
    if isfield(r, 'fista_iteration_limit_reached') ...
            && r.fista_iteration_limit_reached
        fistaText = sprintf('$%d^*$', r.fista_iterations);
    end
    rowArguments = {latexScientific(r.epsilon), ...
        latexScientific(r.sigma), fistaText, r.newton_iterations, ...
        r.pcg_iterations_total, r.total_solver_time, ...
        latexScientific(r.final_pg)};
    if j < numel(records)
        fprintf(fileId, ...
            '$%s$ & $%s$ & %s & %d & %d & %.2f & $%s$ \\\\\n', ...
            rowArguments{:});
    else
        % The enclosing manuscript supplies the final row break after
        % \input, keeping \bottomrule in the surrounding alignment.
        fprintf(fileId, ...
            '$%s$ & $%s$ & %s & %d & %d & %.2f & $%s$\n', ...
            rowArguments{:});
    end
    if j == 3
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

function reportSanity(label, result, targetMass)
massError = abs(src.constraints.Mass( ...
    result.rho, result.grid.h) - targetMass);
fprintf(['  %s sanity: mass %.3e, minrho %.3e, ' ...
    'E %.15e, PG %.3e\n'], label, massError, min(result.rho), ...
    result.target_energy, result.diagnostics.final_pg_residual);
if massError > 1e-12
    warning('%s mass error %.3e exceeds 1e-12.', label, massError);
end
if min(result.rho) < -1e-13
    warning('%s minimum density %.3e violates positivity tolerance.', ...
        label, min(result.rho));
end
end

function reportEnergyIncrease(label, energy)
if numel(energy) < 2
    return;
end
increase = max(diff(energy));
tolerance = 1e-12 * max(1, max(abs(energy)));
if increase > tolerance
    warning('%s has a material recorded energy increase %.3e.', ...
        label, increase);
else
    fprintf('  %s maximum energy increase: %.3e\n', label, increase);
end
end

function text = yesNo(value)
if value
    text = 'yes';
else
    text = 'no';
end
end
