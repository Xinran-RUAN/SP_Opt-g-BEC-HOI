%RUN_REGULARIZATION_COMPARISON Fisher-denominator convergence experiment.
%
% Three convexity-preserving r_epsilon families are compared against one
% shared small-epsilon comparison state. The reported energy error is
% evaluated with the common unregularized Fisher functional r_0(rho)=rho.
clearvars; clc;

% ======================== paper configuration =========================
fisher_choices = {'shift', 'piecewise_c1', 'piecewise_c2'};
epsilon_list = [1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3];
epsilon_reference = 1e-5;

sigma = 1e-2;
beta = 10;
delta = 10;
mass = 1;
L = 32;
N = 512;
R0 = 0.75 * L;
R1 = 0.90 * L;

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-10;
reuse_checkpoint = true;
% =====================================================================

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
output_folder = fullfile(root, 'results', 'fisher_regularization');
figure_folder = fullfile(root, 'results', 'figures');
if ~isfolder(output_folder), mkdir(output_folder); end
if ~isfolder(figure_folder), mkdir(figure_folder); end
checkpoint_file = fullfile(output_folder, ...
    'fisher_regularization_convergence_checkpoint.mat');
data_file = fullfile(output_folder, ...
    'fisher_regularization_convergence_data.mat');
state_fig_file = fullfile(figure_folder, ...
    'kinetic_regularization_state_error.fig');
state_eps_file = fullfile(figure_folder, ...
    'kinetic_regularization_state_error.eps');
energy_fig_file = fullfile(figure_folder, ...
    'kinetic_regularization_energy_error.fig');
energy_eps_file = fullfile(figure_folder, ...
    'kinetic_regularization_energy_error.eps');

signature.fisher_choices = fisher_choices;
signature.epsilon_list = epsilon_list;
signature.epsilon_reference = epsilon_reference;
signature.sigma = sigma;
signature.beta = beta;
signature.delta = delta;
signature.mass = mass;
signature.L = L;
signature.N = N;
signature.R0 = R0;
signature.R1 = R1;
signature.potential = 'harmonic core + C-infinity continuation';
checkpoint = loadCheckpoint(checkpoint_file, signature, ...
    numel(fisher_choices), numel(epsilon_list), reuse_checkpoint);

template = makeTemplate(beta, delta, mass, L, N, R0, R1, sigma, ...
    main_pg_switch_tol, final_pg_tol);
grid = model.SetupGrid1D(template.parameters);
rho_initial = exp(-grid.x .^ 2) / sqrt(pi) ...
    + 0.02 / grid.domain_length;
rho_initial = rho_initial * mass ...
    / src.constraints.Mass(rho_initial, grid.h);

for familyIndex = 1:numel(fisher_choices)
    warmStart = rho_initial;
    for epsilonIndex = 1:numel(epsilon_list)
        epsilon = epsilon_list(epsilonIndex);
        if ~isempty(checkpoint.solutions{familyIndex, epsilonIndex})
            state = checkpoint.solutions{familyIndex, epsilonIndex};
            warmStart = state.rho;
            fprintf('Reused %s, epsilon %.1e (PG %.3e).\n', ...
                fisher_choices{familyIndex}, epsilon, ...
                state.diagnostics.final_pg_residual);
            continue;
        end
        fisher = makeFisher(fisher_choices{familyIndex}, epsilon);
        config = configureCase(template, fisher, epsilon);
        fprintf('Solving %s, epsilon %.1e ...\n', ...
            fisher_choices{familyIndex}, epsilon);
        state = experiments.SolveGroundState(config, warmStart);
        checkpoint.solutions{familyIndex, epsilonIndex} = state;
        warmStart = state.rho;
        saveCheckpoint(checkpoint_file, checkpoint);
    end
end

if isempty(checkpoint.reference)
    % A common, much smaller shift is used only to obtain a stable
    % comparison state. All energy errors below use r_0(rho)=rho.
    fisherReference = makeFisher('shift', epsilon_reference);
    referenceConfig = configureCase( ...
        template, fisherReference, epsilon_reference);
    referenceStart = checkpoint.solutions{1, end}.rho;
    fprintf('Solving shared reference, epsilon_ref %.1e ...\n', ...
        epsilon_reference);
    checkpoint.reference = experiments.SolveGroundState( ...
        referenceConfig, referenceStart);
    saveCheckpoint(checkpoint_file, checkpoint);
else
    fprintf('Reused shared epsilon_ref %.1e comparison state.\n', ...
        epsilon_reference);
end
reference = checkpoint.reference;

numberOfFamilies = numel(fisher_choices);
numberOfEpsilons = numel(epsilon_list);
state_error = nan(numberOfFamilies, numberOfEpsilons);
energy_error = nan(numberOfFamilies, numberOfEpsilons);
final_pg = nan(numberOfFamilies, numberOfEpsilons);
min_rho = nan(numberOfFamilies, numberOfEpsilons);
target_energy = nan(numberOfFamilies, numberOfEpsilons);
rho = cell(numberOfFamilies, numberOfEpsilons);
reference_unregularized_energy = unregularizedFisherEnergy( ...
    reference.rho, reference.problem);
for familyIndex = 1:numberOfFamilies
    for epsilonIndex = 1:numberOfEpsilons
        state = checkpoint.solutions{familyIndex, epsilonIndex};
        rho{familyIndex, epsilonIndex} = state.rho;
        difference = state.rho - reference.rho;
        state_error(familyIndex, epsilonIndex) = sqrt( ...
            grid.h * sum(difference .^ 2));
        commonEnergy = unregularizedFisherEnergy( ...
            state.rho, reference.problem);
        energy_error(familyIndex, epsilonIndex) = abs( ...
            commonEnergy - reference_unregularized_energy);
        final_pg(familyIndex, epsilonIndex) = ...
            state.diagnostics.final_pg_residual;
        min_rho(familyIndex, epsilonIndex) = min(state.rho);
        target_energy(familyIndex, epsilonIndex) = state.target_energy;
    end
end
state_rate = epsilonRates(epsilon_list, state_error);
energy_rate = epsilonRates(epsilon_list, energy_error);

labels = {'Shift', 'Piecewise $C^1$', 'Piecewise $C^2$'};
[state_reference_lines] = makeFigure( ...
    epsilon_list, state_error, labels, ...
    'state error', 'Kinetic-regularization state error', ...
    state_fig_file, state_eps_file);
[energy_reference_lines] = makeFigure( ...
    epsilon_list, energy_error, labels, ...
    'energy error', 'Kinetic-regularization energy error', ...
    energy_fig_file, energy_eps_file);
printTables(fisher_choices, epsilon_list, state_error, energy_error, ...
    state_rate, energy_rate, final_pg, min_rho);

metadata = signature;
metadata.reference_description = ...
    'shared shift state at epsilon_ref; errors use r_0(rho)=rho';
metadata.state_error_definition = ...
    'sqrt(h*sum((rho_epsilon-rho_ref).^2))';
metadata.energy_error_definition = ...
    '|E_ref(rho_epsilon)-E_ref(rho_ref)| with r_0(rho)=rho';
metadata.reference_final_pg = reference.diagnostics.final_pg_residual;
metadata.reference_min_rho = min(reference.rho);
reference_rho = reference.rho;
save(data_file, 'metadata', 'fisher_choices', 'labels', ...
    'epsilon_list', 'epsilon_reference', 'sigma', 'beta', 'delta', ...
    'mass', 'L', 'N', 'R0', 'R1', 'rho', 'reference_rho', ...
    'reference_unregularized_energy', 'state_error', 'energy_error', ...
    'state_rate', 'energy_rate', 'final_pg', 'min_rho', ...
    'target_energy', 'state_reference_lines', ...
    'energy_reference_lines', '-v7.3');

fprintf('\nSaved Fisher-regularization outputs\n');
fprintf('  %s\n', state_fig_file);
fprintf('  %s\n', state_eps_file);
fprintf('  %s\n', energy_fig_file);
fprintf('  %s\n', energy_eps_file);
fprintf('  %s\n', data_file);

function template = makeTemplate(beta, delta, mass, L, N, R0, R1, ...
    sigma, mainTolerance, finalTolerance)
template.parameters = model.DefaultParameters1D();
template.parameters.beta = beta;
template.parameters.delta = delta;
template.parameters.mass = mass;
template.parameters.L = L;
template.parameters.N = N;
template.R0 = R0;
template.R1 = R1;
template.sigma = sigma;
template.main_tolerance = mainTolerance;
template.final_tolerance = finalTolerance;
end

function config = configureCase(template, fisher, epsilon)
config = experiments.DefaultConfig();
config.parameters = template.parameters;
config.parameters.epsilon = epsilon;
config.fisher_regularization = fisher;
config.trapping_potential.V = @(x) ...
    model.HarmonicCInfPeriodicPotential( ...
        x, template.parameters.L, template.R0, template.R1);
config.trapping_potential.label = ...
    'harmonic core + C-infinity far-field continuation';
sigma = template.sigma;
config.potential_regularization.name = 'inline_p_sigma';
config.potential_regularization.label = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma';
config.potential_regularization.sigma = sigma;
config.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
config.potential_regularization.dp_sigma = @(rho) ...
    rho ./ hypot(rho, sigma);
config.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
config.potential_regularization.prox_type = 'generic_convex';
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = 'semismooth';
config.solver.projection_tol = 1e-13;
config.solver.pg_tol = template.main_tolerance;
config.solver.final_pg_tol = template.final_tolerance;
config.solver.certification_tol = template.final_tolerance;
config.solver.max_iter = 200000;
config.solver.display = false;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = template.main_tolerance;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = template.final_tolerance;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.potential_prox.mass_tol = 1e-12;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.potential_prox.inner_max_iter = 80;
end

function fisher = makeFisher(choice, epsilon)
switch lower(choice)
    case 'shift'
        % Simple formulas remain visible and directly editable here.
        fisher.r_epsilon = @(rho) rho + epsilon;
        fisher.dr_epsilon = @(rho) ones(size(rho));
        fisher.d2r_epsilon = @(rho) zeros(size(rho));
        fisher.epsilon = epsilon;
        fisher.label = 'r_epsilon(rho)=rho+epsilon';
        fisher.name = 'shift';
        fisher.regularity = inf;
    case 'piecewise_c1'
        fisher = src.regularization.PiecewiseFisherCm(epsilon, 1);
    case 'piecewise_c2'
        fisher = src.regularization.PiecewiseFisherCm(epsilon, 2);
    otherwise
        error('Unknown Fisher choice "%s".', choice);
end
fisher = src.regularization.NormalizeFisher(fisher);
end

function energy = unregularizedFisherEnergy(rho, problem)
rho = rho(:);
if any(rho <= 0)
    error('Reference Fisher energy requires a strictly positive state.');
end
q = src.discretization.ps.FirstDerivative(rho, problem.plan);
potential = problem.potential_regularization;
p = potential.p_sigma(rho);
energy = problem.grid.h * sum(q .^ 2 ./ (8 * rho) ...
    + problem.V(:) .* p + problem.beta / 2 * rho .^ 2 ...
    + problem.delta / 2 * q .^ 2);
end

function rates = epsilonRates(epsilon, errors)
rates = nan(size(errors));
for row = 1:size(errors, 1)
    for column = 2:size(errors, 2)
        if errors(row, column - 1) > 1e-14 ...
                && errors(row, column) > 1e-14
            rates(row, column) = log(errors(row, column - 1) ...
                / errors(row, column)) / log(epsilon(column - 1) ...
                / epsilon(column));
        end
    end
end
end

function referenceLines = makeFigure( ...
    epsilon, errors, labels, yLabel, plotTitle, ...
    figFile, epsFile)
fontSize = 20;
lineWidth = 2;
markerSize = 9;
styles = {'-o', '--s', '-.^'};
colors = lines(3);
f = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 650]);
ax = axes('Parent', f);
hold(ax, 'on');
for row = 1:size(errors, 1)
    loglog(ax, epsilon, errors(row, :), styles{row}, ...
        'Color', colors(row, :), 'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, 'DisplayName', labels{row});
end

% Plot first- and second-order slope guides.  Their vertical offsets are
% chosen from the data scale at epsilon=1e-2; only the slopes are meant
% to be interpreted.
[~, anchorIndex] = min(abs(log(epsilon) - log(1e-2)));
anchorEpsilon = epsilon(anchorIndex);
anchorError = median(errors(:, anchorIndex));
referenceLines.epsilon = epsilon;
referenceLines.anchor_epsilon = anchorEpsilon;
referenceLines.first_order = 2 * anchorError ...
    * (epsilon / anchorEpsilon);
referenceLines.second_order = 0.5 * anchorError ...
    * (epsilon / anchorEpsilon) .^ 2;
loglog(ax, epsilon, referenceLines.first_order, '--', ...
    'Color', [0.10, 0.10, 0.10], 'LineWidth', lineWidth, ...
    'DisplayName', '$\mathcal{O}(\varepsilon)$');
loglog(ax, epsilon, referenceLines.second_order, ':', ...
    'Color', [0.50, 0.50, 0.50], 'LineWidth', lineWidth, ...
    'DisplayName', '$\mathcal{O}(\varepsilon^2)$');
set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse', ...
    'FontSize', fontSize, 'LineWidth', 1);
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, '$\varepsilon$', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
ylabel(ax, yLabel, 'Interpreter', 'latex', 'FontSize', fontSize);
title(ax, plotTitle, 'FontSize', fontSize, 'FontWeight', 'normal');
legend(ax, 'Location', 'best', 'Interpreter', 'latex', ...
    'FontSize', fontSize);
set(f, 'PaperPositionMode', 'auto', 'Visible', 'on');
savefig(f, figFile);
set(f, 'Visible', 'off', 'Renderer', 'painters');
print(f, epsFile, '-depsc2', '-vector');
close(f);
end

function printTables(names, epsilon, stateError, energyError, ...
    stateRate, energyRate, finalPg, minRho)
fprintf('\n============================================================\n');
fprintf('Fisher regularization convergence\n');
fprintf('============================================================\n');
for row = 1:numel(names)
    fprintf('\n%s\n', names{row});
    fprintf(' epsilon    state_error  rate    energy_error rate    PG       minrho\n');
    for column = 1:numel(epsilon)
        fprintf(' %.1e  %.6e  %5s  %.6e  %5s  %.2e  %.3e\n', ...
            epsilon(column), stateError(row, column), ...
            rateText(stateRate(row, column)), energyError(row, column), ...
            rateText(energyRate(row, column)), finalPg(row, column), ...
            minRho(row, column));
    end
end
end

function text = rateText(value)
if isnan(value)
    text = '--';
else
    text = sprintf('%.2f', value);
end
end

function checkpoint = loadCheckpoint(file, signature, rows, columns, reuse)
checkpoint.signature = signature;
checkpoint.solutions = cell(rows, columns);
checkpoint.reference = [];
if ~reuse || ~isfile(file)
    return;
end
loaded = load(file, 'checkpoint');
if isfield(loaded, 'checkpoint') ...
        && isequaln(loaded.checkpoint.signature, signature) ...
        && isequal(size(loaded.checkpoint.solutions), [rows, columns])
    checkpoint = loaded.checkpoint;
    fprintf('Reusing compatible Fisher convergence checkpoint.\n');
else
    fprintf('Ignoring incompatible Fisher convergence checkpoint.\n');
end
end

function saveCheckpoint(file, checkpoint)
save(file, 'checkpoint', '-v7.3');
end
