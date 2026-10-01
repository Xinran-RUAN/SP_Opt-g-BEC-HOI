%RUN_PLOT_SIGMA_SPECTRAL_ACCURACY Effect of sigma on Fourier accuracy.
%
% Each sigma defines one fixed regularized problem.  The script compares
% mesh-refinement errors only when epsilon, beta, delta, L, N_list, N_ref,
% and the trapping potential agree.  State and energy errors are exported
% as separate paper figures (no subplot/tiledlayout).
clearvars; clc;

% ======================== paper configuration =========================
epsilon = 1e-2;
sigma_list = [0, 1e-2, 1e-4, 1e-6];
beta = 10;
delta = 10;
mass = 1;
L = 32;
N_list = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, ...
    16384, 32768];
N_ref = 65536;

R0_fraction = 0.75;
R1_fraction = 0.90;
R0 = R0_fraction * L;
R1 = R1_fraction * L;
potential_choice = 'harmonic_cinf_periodic';
potential_label = [ ...
    'harmonic core plus far-field C-infinity periodic continuation'];

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
reuse_results = true;
solve_missing_cases = true;
% =====================================================================

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();

output_folder = fullfile(root, 'results', 'figures');
case_folder = fullfile(root, 'results', 'sigma_spectral_accuracy');
paper_figure_folder = fullfile(root, 'figs');
repository_root = fileparts(root);
manuscript_figure_folder = fullfile(repository_root, 'manuscript', 'figs');
if ~isfolder(output_folder), mkdir(output_folder); end
if ~isfolder(case_folder), mkdir(case_folder); end
if ~isfolder(paper_figure_folder), mkdir(paper_figure_folder); end

state_figure_file = fullfile(output_folder, ...
    'sigma_effect_state_error_1d.fig');
state_eps_file = fullfile(output_folder, ...
    'sigma_effect_state_error_1d.eps');
energy_figure_file = fullfile(output_folder, ...
    'sigma_effect_energy_error_1d.fig');
energy_eps_file = fullfile(output_folder, ...
    'sigma_effect_energy_error_1d.eps');
data_file = fullfile(root, 'results', 'sigma_effect_with_sigma0.mat');

validateMeshes(N_list, N_ref);
solve_N = [N_list, N_ref];
validateMeshes(solve_N(1:end-1), N_ref);

template.epsilon = epsilon;
template.beta = beta;
template.delta = delta;
template.mass = mass;
template.L = L;
template.R0 = R0;
template.R1 = R1;
template.R0_fraction = R0_fraction;
template.R1_fraction = R1_fraction;
template.potential_choice = potential_choice;
template.potential_label = potential_label;
template.main_pg_switch_tol = main_pg_switch_tol;
template.final_pg_tol = final_pg_tol;
template.s_epsilon = @(rho) rho + epsilon;
template.ds_epsilon = @(rho) ones(size(rho));
template.d2s_epsilon = @(rho) zeros(size(rho));
template.fisher_label = 's_epsilon(rho)=rho+epsilon';
template.V = @(x) model.HarmonicCInfPeriodicPotential(x, L, R0, R1);

% These are only preferred data sources.  Every file is checked before
% use; the pre-existing sigma=1e-4 file is deliberately listed so that an
% incompatible periodic-cosine trap is diagnosed and rejected.
preferred_files = cell(size(sigma_list));
preferred_files{1} = '';
preferred_files{2} = fullfile(root, 'results', 'potential_mesh', ...
    'potential_mesh_inline_p_sigma_sigma_1em02_Nref8192_eps_1em02.mat');
preferred_files{3} = fullfile(root, 'results', 'potential_mesh', ...
    'potential_mesh_inline_p_sigma_sigma_1em04_Nref8192_eps_1em02.mat');
preferred_files{4} = '';

case_data = cell(numel(sigma_list), 1);
for sigma_index = 1:numel(sigma_list)
    sigma = sigma_list(sigma_index);
    case_file = fullfile(case_folder, sprintf( ...
        'sigma_mesh_%s_eps1em02_L32_Nref%d.mat', ...
        numberTag(sigma), N_ref));
    checkpoint_file = fullfile(case_folder, sprintf( ...
        'sigma_mesh_%s_eps1em02_L32_checkpoint.mat', ...
        numberTag(sigma)));

    loaded = false;
    if reuse_results && isfile(case_file)
        [loaded, candidate, reason] = loadCompatibleCase( ...
            case_file, sigma, template, N_list, N_ref);
        if loaded
            case_data{sigma_index} = candidate;
            fprintf('sigma %.1e: reused %s\n', sigma, case_file);
        else
            fprintf('sigma %.1e: rejected %s (%s).\n', ...
                sigma, case_file, reason);
        end
    end
    if ~loaded && reuse_results && ~isempty(preferred_files{sigma_index}) ...
            && isfile(preferred_files{sigma_index})
        [loaded, candidate, reason] = loadCompatibleCase( ...
            preferred_files{sigma_index}, sigma, template, N_list, N_ref);
        if loaded
            case_data{sigma_index} = candidate;
            fprintf('sigma %.1e: reused compatible existing mesh data.\n', ...
                sigma);
        else
            fprintf(['sigma %.1e: existing MAT is incompatible and was ' ...
                'not mixed into the plot (%s).\n'], sigma, reason);
        end
    end

    if ~loaded
        if ~solve_missing_cases
            error('No compatible fixed-sigma mesh data for sigma %.3e.', ...
                sigma);
        end
        fprintf('sigma %.1e: solving missing compatible mesh case.\n', sigma);
        seed_solutions = loadSeedSolutions( ...
            root, sigma, solve_N, epsilon, L);
        [solutions, solved_case] = solveSigmaCase(sigma, template, ...
            solve_N, N_list, N_ref, checkpoint_file, seed_solutions);
        case_data{sigma_index} = solved_case;
        caseData = solved_case; %#ok<NASGU>
        save(case_file, 'caseData', '-v7.3');
    end
end

state_error_curves = cell2mat(cellfun(@(entry) entry.state_error, ...
    case_data, 'UniformOutput', false));
energy_error_curves = cell2mat(cellfun(@(entry) entry.energy_error, ...
    case_data, 'UniformOutput', false));
final_pg_curves = cell2mat(cellfun(@(entry) entry.final_pg, ...
    case_data, 'UniformOutput', false));
reference_pg = cellfun(@(entry) entry.reference_pg, case_data);
state_rate_curves = adjacentRates(state_error_curves);
energy_rate_curves = adjacentRates(energy_error_curves);
labels = arrayfun(@sigmaLabel, sigma_list, 'UniformOutput', false);
source_files = cellfun(@(entry) entry.source_file, case_data, ...
    'UniformOutput', false);
certified = cellfun(@(entry) entry.certified, case_data);

metadata.epsilon = epsilon;
metadata.beta = beta;
metadata.delta = delta;
metadata.mass = mass;
metadata.L = L;
metadata.N_ref = N_ref;
metadata.potential_choice = potential_choice;
metadata.potential_label = potential_label;
metadata.R0_fraction = R0_fraction;
metadata.R1_fraction = R1_fraction;
metadata.R0 = R0;
metadata.R1 = R1;
metadata.state_error_definition = 'total_spectral_L2_error';
metadata.energy_error_definition = ...
    '|E_N(rho_N)-E_Nref(rho_ref)|';

save(data_file, 'sigma_list', 'N_list', 'N_ref', ...
    'state_error_curves', 'energy_error_curves', 'final_pg_curves', ...
    'state_rate_curves', 'energy_rate_curves', 'reference_pg', ...
    'labels', 'metadata', 'source_files', 'certified', 'case_data', ...
    '-v7.3');

makeSeparateFigures(N_list, sigma_list, state_error_curves, ...
    energy_error_curves, state_figure_file, state_eps_file, ...
    energy_figure_file, energy_eps_file);
syncFigureFiles({state_figure_file, state_eps_file, ...
    replace(state_eps_file, '.eps', '.png'), energy_figure_file, ...
    energy_eps_file, replace(energy_eps_file, '.eps', '.png')}, ...
    paper_figure_folder);
if isfolder(manuscript_figure_folder)
    syncFigureFiles({state_eps_file, ...
        replace(state_eps_file, '.eps', '.png'), energy_eps_file, ...
        replace(energy_eps_file, '.eps', '.png')}, ...
        manuscript_figure_folder);
end

fprintf('\nSigma effect on Fourier spatial accuracy\n');
fprintf('  epsilon=%g, beta=%g, delta=%g, L=%g, Nref=%d\n', ...
    epsilon, beta, delta, L, N_ref);
for sigma_index = 1:numel(sigma_list)
    fprintf('  sigma=%-9.1e N =', sigma_list(sigma_index));
    fprintf(' %d', N_list);
    fprintf('   certified=%d\n', certified(sigma_index));
end
fprintf('State figure : %s\n', state_figure_file);
fprintf('State EPS    : %s\n', state_eps_file);
fprintf('Energy figure: %s\n', energy_figure_file);
fprintf('Energy EPS   : %s\n', energy_eps_file);
fprintf('Data MAT     : %s\n', data_file);

function [solutions, caseData] = solveSigmaCase(sigma, template, ...
    solveN, reportN, Nref, checkpointFile, seedSolutions)
signature = [template.epsilon, sigma, template.beta, template.delta, ...
    template.mass, template.L, template.R0, template.R1, solveN];
solutions = cell(numel(solveN), 1);
if isfile(checkpointFile)
    saved = load(checkpointFile, 'checkpoint');
    if isfield(saved, 'checkpoint') ...
            && isfield(saved.checkpoint, 'signature') ...
            && isequal(saved.checkpoint.signature, signature)
        solutions = saved.checkpoint.solutions;
        fprintf('  reused compatible checkpoint.\n');
    elseif isfield(saved, 'checkpoint') ...
            && isfield(saved.checkpoint, 'signature') ...
            && numel(saved.checkpoint.signature) >= 9 ...
            && isequal(saved.checkpoint.signature(1:8), signature(1:8))
        oldN = saved.checkpoint.signature(9:end);
        solutions = mergeSolutions(solutions, solveN, ...
            saved.checkpoint.solutions, oldN);
        fprintf('  reused compatible grids from an older checkpoint.\n');
    else
        fprintf('  ignored incompatible checkpoint.\n');
    end
end
solutions = mergeSolutions(solutions, solveN, seedSolutions, solveN);

previous = [];
for index = 1:numel(solveN)
    N = solveN(index);
    if ~isempty(solutions{index})
        previous = solutions{index};
        fprintf('  N=%d reused (PG %.3e).\n', N, ...
            previous.diagnostics.final_pg_residual);
        continue;
    end
    config = makeConfig(N, sigma, template);
    grid = model.SetupGrid1D(config.parameters);
    if isempty(previous)
        rho0 = seedAtN(seedSolutions, solveN, N);
        if isempty(rho0)
            rho0 = exp(-grid.x .^ 2);
            rho0 = rho0 * (template.mass ...
                / src.constraints.Mass(rho0, grid.h));
        end
        rho0 = src.constraints.ProjectPositiveConservative( ...
            rho0, template.mass, grid.h, config.solver.projection_tol);
    else
        transferProblem = previous.problem;
        transferProblem.grid = grid;
        transferProblem.plan = src.discretization.ps.Plan1D(grid);
        transferProblem.V = template.V(grid.x);
        transferSolver = config.solver;
        transferSolver.splitting = 'legacy_full_gradient';
        rho0 = experiments.TransferState1D(previous.rho, grid, ...
            template.mass, transferProblem, transferSolver);
    end
    timer = tic;
    result = experiments.SolveGroundState(config, rho0);
    result.solve_wall_time = toc(timer);
    solutions{index} = result;
    previous = result;
    fprintf('  N=%d done %.2fs, E %.15e, PG %.3e\n', ...
        N, result.solve_wall_time, result.target_energy, ...
        result.diagnostics.final_pg_residual);
    checkpoint.signature = signature;
    checkpoint.solutions = solutions;
    save(checkpointFile, 'checkpoint', '-v7.3');
end

reference = solutionAt(solutions, solveN, Nref);
stateError = zeros(1, numel(reportN));
energyError = zeros(1, numel(reportN));
finalPG = zeros(1, numel(reportN));
for index = 1:numel(reportN)
    state = solutionAt(solutions, solveN, reportN(index));
    comparison = src.diagnostics.SpectralStateComparison( ...
        state.rho, state.grid, reference.rho, reference.grid);
    stateError(index) = comparison.total_spectral_L2_error;
    energyError(index) = abs(state.target_energy - reference.target_energy);
    finalPG(index) = state.diagnostics.final_pg_residual;
end
caseData = buildCaseData(sigma, template, reportN, Nref, ...
    stateError, energyError, finalPG, checkpointFile, ...
    reference.diagnostics.final_pg_residual);
end

function seeds = loadSeedSolutions(root, sigma, solveN, epsilon, L)
seeds = cell(numel(solveN), 1);
if abs(epsilon - 1e-2) > eps || L ~= 32
    return;
end
if abs(sigma - 1e-2) <= eps
    filename = fullfile(root, 'results', 'figures', ...
        'mesh_convergence_1d_checkpoint.mat');
    if ~isfile(filename), return; end
    loaded = load(filename, 'checkpoint');
    if isfield(loaded.checkpoint, 'continued') ...
            && isfield(loaded.checkpoint.continued, 'N_values') ...
            && isfield(loaded.checkpoint.continued, 'solutions')
        seeds = mergeSolutions(seeds, solveN, ...
            loaded.checkpoint.continued.solutions, ...
            loaded.checkpoint.continued.N_values);
    end
    return;
end
legacy = fullfile(root, 'results', 'sigma_spectral_accuracy', ...
    sprintf(['sigma_mesh_%s_eps1em02_L32_Nref8192_' ...
    'checkpoint.mat'], numberTag(sigma)));
if ~isfile(legacy), return; end
loaded = load(legacy, 'checkpoint');
if ~isfield(loaded.checkpoint, 'signature') ...
        || numel(loaded.checkpoint.signature) < 9 ...
        || ~isfield(loaded.checkpoint, 'solutions')
    return;
end
oldN = loaded.checkpoint.signature(9:end);
seeds = mergeSolutions(seeds, solveN, ...
    loaded.checkpoint.solutions, oldN);
end

function destination = mergeSolutions(destination, destinationN, source, sourceN)
if isempty(source) || ~iscell(source)
    return;
end
for index = 1:numel(destinationN)
    sourceIndex = find(sourceN == destinationN(index), 1);
    if ~isempty(sourceIndex) && sourceIndex <= numel(source) ...
            && ~isempty(source{sourceIndex}) && isempty(destination{index})
        destination{index} = source{sourceIndex};
    end
end
end

function config = makeConfig(N, sigma, template)
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
if sigma == 0
    % Use the exact unsmoothed linear potential contribution.  In
    % particular, do not evaluate the square-root formula at sigma=0,
    % where its derivative representation is ambiguous at vacuum.
    config.potential_regularization = src.potential.MakeLinear();
else
    config.potential_regularization.name = 'inline_fixed_sigma';
    config.potential_regularization.sigma = sigma;
    config.potential_regularization.p_sigma = @(rho) ...
        rho .^ 2 ./ (hypot(rho, sigma) + sigma);
    config.potential_regularization.dp_sigma = @(rho) ...
        rho ./ hypot(rho, sigma);
    config.potential_regularization.d2p_sigma = @(rho) ...
        (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
    config.potential_regularization.label = sprintf( ...
        'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma, sigma=%.3e', sigma);
    config.potential_regularization.prox_type = 'generic_convex';
    config.potential_regularization.convexity_tol = 1e-14;
end
config.trapping_potential.V = template.V;
config.trapping_potential.label = template.potential_label;
config.trapping_potential.mode = template.potential_choice;
config.trapping_potential.boundary_periodicized = true;
config.trapping_potential.modification_start = template.R0;
config.trapping_potential.transition_end = template.R1;
config.trapping_potential.reference_V = @(x) 0.5 * x .^ 2;
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

function [ok, caseData, reason] = loadCompatibleCase(filename, sigma, ...
    template, reportN, Nref)
ok = false;
caseData = struct();
reason = '';
loaded = load(filename);
if isfield(loaded, 'caseData')
    candidate = loaded.caseData;
    required = {'sigma', 'epsilon', 'beta', 'delta', 'mass', 'L', ...
        'N_list', 'N_ref', 'potential_choice', 'state_error', ...
        'energy_error', 'final_pg'};
    if ~all(isfield(candidate, required))
        reason = 'dedicated case file is incomplete';
        return;
    end
    if ~sameScalar(candidate.sigma, sigma) ...
            || ~sameScalar(candidate.epsilon, template.epsilon) ...
            || ~sameScalar(candidate.beta, template.beta) ...
            || ~sameScalar(candidate.delta, template.delta) ...
            || ~sameScalar(candidate.mass, template.mass) ...
            || ~sameScalar(candidate.L, template.L) ...
            || candidate.N_ref ~= Nref ...
            || ~isequal(candidate.N_list(:).', reportN) ...
            || ~strcmpi(candidate.potential_choice, ...
                template.potential_choice)
        reason = 'dedicated metadata mismatch';
        return;
    end
    caseData = candidate;
    caseData.source_file = filename;
    ok = true;
    return;
end

requiredTop = {'parameters', 'potential_regularization', ...
    'trapping_potential', 'diagnostics'};
if ~all(isfield(loaded, requiredTop)) ...
        || ~isfield(loaded.diagnostics, 'records')
    reason = 'archive fields missing';
    return;
end
parameters = loaded.parameters;
potential = loaded.potential_regularization;
trap = loaded.trapping_potential;
if ~isfield(parameters, 'epsilon') || ~isfield(parameters, 'L') ...
        || ~isfield(parameters, 'beta') || ~isfield(parameters, 'delta') ...
        || ~isfield(parameters, 'mass') || ~isfield(parameters, 'N_ref') ...
        || ~isfield(potential, 'sigma')
    reason = 'essential archive metadata missing';
    return;
end
if ~sameScalar(parameters.epsilon, template.epsilon) ...
        || ~sameScalar(parameters.L, template.L) ...
        || ~sameScalar(parameters.beta, template.beta) ...
        || ~sameScalar(parameters.delta, template.delta) ...
        || ~sameScalar(parameters.mass, template.mass) ...
        || parameters.N_ref ~= Nref || ~sameScalar(potential.sigma, sigma)
    reason = 'numeric metadata mismatch';
    return;
end
if ~isfield(trap, 'mode') ...
        || ~strcmpi(char(trap.mode), template.potential_choice)
    reason = 'trapping potential mismatch';
    return;
end
records = loaded.diagnostics.records;
recordN = [records.N];
if ~all(ismember(reportN, recordN))
    reason = 'required mesh records missing';
    return;
end
stateError = zeros(1, numel(reportN));
energyError = zeros(1, numel(reportN));
finalPG = zeros(1, numel(reportN));
for index = 1:numel(reportN)
    record = records(find(recordN == reportN(index), 1));
    stateError(index) = record.total_spectral_L2_error;
    energyError(index) = record.target_energy_error;
    finalPG(index) = record.final_pg_residual;
end
referencePG = NaN;
if isfield(loaded.diagnostics, 'finest') ...
        && isfield(loaded.diagnostics.finest, 'final_pg_residual')
    referencePG = loaded.diagnostics.finest.final_pg_residual;
end
caseData = buildCaseData(sigma, template, reportN, Nref, ...
    stateError, energyError, finalPG, filename, referencePG);
ok = true;
end

function caseData = buildCaseData(sigma, template, Nlist, Nref, ...
    stateError, energyError, finalPG, sourceFile, referencePG)
caseData.sigma = sigma;
caseData.epsilon = template.epsilon;
caseData.beta = template.beta;
caseData.delta = template.delta;
caseData.mass = template.mass;
caseData.L = template.L;
caseData.N_list = Nlist;
caseData.N_ref = Nref;
caseData.potential_choice = template.potential_choice;
caseData.potential_label = template.potential_label;
caseData.R0 = template.R0;
caseData.R1 = template.R1;
caseData.state_error = stateError;
caseData.energy_error = energyError;
caseData.final_pg = finalPG;
caseData.reference_pg = referencePG;
caseData.certified = all(isfinite(finalPG)) ...
    && all(finalPG <= 1e-9) && isfinite(referencePG) ...
    && referencePG <= 1e-9;
caseData.source_file = sourceFile;
end

function rho = seedAtN(seedSolutions, solveN, N)
rho = [];
if isempty(seedSolutions) || ~iscell(seedSolutions)
    return;
end
index = find(solveN == N, 1);
if isempty(index) || index > numel(seedSolutions) ...
        || isempty(seedSolutions{index})
    return;
end
candidate = seedSolutions{index};
if isstruct(candidate) && isfield(candidate, 'rho')
    rho = candidate.rho;
end
end

function state = solutionAt(solutions, Nvalues, N)
index = find(Nvalues == N, 1);
if isempty(index) || isempty(solutions{index})
    error('Missing computed state at N=%d.', N);
end
state = solutions{index};
end

function makeSeparateFigures(Nvalues, sigmaValues, stateErrors, ...
    energyErrors, stateFigFile, stateEpsFile, energyFigFile, ...
    energyEpsFile)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
lineStyles = {'-o', '--s', '-.^', ':d', '-x', '--v', '-.p', ':h'};
colors = lines(numel(sigmaValues));

stateFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 650], ...
    'Name', 'Effect of sigma on Fourier state accuracy');
stateAxis = axes(stateFigure);
set(stateAxis, 'XScale', 'log', 'YScale', 'log');
hold(stateAxis, 'on');
for index = 1:numel(sigmaValues)
    loglog(stateAxis, Nvalues, stateErrors(index, :), ...
        lineStyles{1 + mod(index - 1, numel(lineStyles))}, ...
        'Color', colors(index, :), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'DisplayName', sigmaLabel(sigmaValues(index)));
end
grid(stateAxis, 'on'); box(stateAxis, 'on');
xlabel(stateAxis, '$N$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(stateAxis, '$e_{\rho,2}$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
title(stateAxis, 'Effect of $\sigma$ on Fourier state accuracy', ...
    'Interpreter', 'latex', 'FontSize', fontSize, ...
    'FontWeight', 'normal');
legend(stateAxis, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(stateAxis, 'FontSize', fontSize, 'LineWidth', 1, ...
    'XTick', Nvalues, 'XLim', [0.85 * min(Nvalues), 1.18 * max(Nvalues)]);
set(stateFigure, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(stateFigure, stateFigFile);
set(stateFigure, 'Visible', 'off');
print(stateFigure, stateEpsFile, '-depsc2', '-vector');
exportgraphics(stateFigure, replace(stateEpsFile, '.eps', '.png'), ...
    'Resolution', 300);
close(stateFigure);

energyFigure = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 650], ...
    'Name', 'Effect of sigma on Fourier energy accuracy');
energyAxis = axes(energyFigure);
set(energyAxis, 'XScale', 'log', 'YScale', 'log');
hold(energyAxis, 'on');
for index = 1:numel(sigmaValues)
    loglog(energyAxis, Nvalues, energyErrors(index, :), ...
        lineStyles{1 + mod(index - 1, numel(lineStyles))}, ...
        'Color', colors(index, :), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'DisplayName', sigmaLabel(sigmaValues(index)));
end
grid(energyAxis, 'on'); box(energyAxis, 'on');
xlabel(energyAxis, '$N$', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel(energyAxis, '$e_E$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
title(energyAxis, 'Effect of $\sigma$ on Fourier energy accuracy', ...
    'Interpreter', 'latex', 'FontSize', fontSize, ...
    'FontWeight', 'normal');
legend(energyAxis, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(energyAxis, 'FontSize', fontSize, 'LineWidth', 1, ...
    'XTick', Nvalues, 'XLim', [0.85 * min(Nvalues), 1.18 * max(Nvalues)]);
set(energyFigure, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(energyFigure, energyFigFile);
set(energyFigure, 'Visible', 'off');
print(energyFigure, energyEpsFile, '-depsc2', '-vector');
exportgraphics(energyFigure, replace(energyEpsFile, '.eps', '.png'), ...
    'Resolution', 300);
close(energyFigure);
end

function label = sigmaLabel(sigma)
if sigma == 0
    label = '$\sigma=0$';
    return;
end
exponent = round(log10(sigma));
if abs(sigma - 10 ^ exponent) <= 100 * eps(max(1, sigma))
    label = sprintf('$\\sigma=10^{%d}$', exponent);
else
    label = sprintf('$\\sigma=%.3g$', sigma);
end
end

function rates = adjacentRates(errors)
rates = NaN(size(errors));
rates(:, 2:end) = log2(errors(:, 1:end-1) ./ errors(:, 2:end));
end

function validateMeshes(Nvalues, Nreference)
if any(mod(Nvalues, 2) ~= 0) || mod(Nreference, 2) ~= 0 ...
        || Nreference <= max(Nvalues) ...
        || any(mod(Nreference, Nvalues) ~= 0) ...
        || any(diff(Nvalues) <= 0)
    error('Use strictly increasing nested even Fourier grids.');
end
end

function tf = sameScalar(a, b)
tf = isscalar(a) && isscalar(b) && isfinite(a) && isfinite(b) ...
    && abs(a - b) <= 100 * eps(max([1, abs(a), abs(b)]));
end

function syncFigureFiles(files, destination)
for index = 1:numel(files)
    [~, name, extension] = fileparts(files{index});
    copyfile(files{index}, fullfile(destination, [name, extension]));
end
end

function tag = numberTag(value)
if value == 0
    tag = '0';
    return;
end
tag = strrep(strrep(sprintf('%.0e', value), '-', 'm'), '+', 'p');
end
