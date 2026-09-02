%RUN_DIAG_HARMONIC_BOUNDARY_SENSITIVITY Fixed-h box sensitivity study.
%
% Direct harmonic states are sensitivity references only. This script does
% not perform a mesh-refinement study and changes no production solver.
clearvars; clc;

epsilon = 1e-2;
sigma = 2e-4;
beta = 10;
delta = 10;
mass = 1;
L_list = [8, 10, 12];
L_base = 8;
N_base = 1024;
h_target = 2 * L_base / N_base;
R_core = 5;
target_state_accuracy = 1e-8;
target_energy_accuracy = 1e-10;

s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho)=rho+epsilon';
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
potential_label = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma, fixed sigma=2e-4';

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
outputFolder = fullfile(root, 'results', ...
    'harmonic_boundary_sensitivity');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputFile = fullfile(outputFolder, ...
    'harmonic_boundary_sensitivity_eps1em2_sigma2em4.mat');

solver = makeSolver();
template.epsilon = epsilon;
template.sigma = sigma;
template.beta = beta;
template.delta = delta;
template.mass = mass;
template.s_epsilon = s_epsilon;
template.ds_epsilon = ds_epsilon;
template.d2s_epsilon = d2s_epsilon;
template.fisher_label = fisher_label;
template.p_sigma = p_sigma;
template.dp_sigma = dp_sigma;
template.d2p_sigma = d2p_sigma;
template.potential_label = potential_label;

caseCount = numel(L_list);
records = struct([]);
harmonicResults = cell(caseCount, 1);
periodicResults = cell(caseCount, 1);
potentialShapes = cell(caseCount, 1);

for index = 1:caseCount
    record = struct();
    L = L_list(index);
    N = 2 * round((2 * L / h_target) / 2);
    grid = makeGrid(L, N);
    R0 = 0.75 * L;
    R1 = 0.90 * L;
    H = 0.5 * grid.x .^ 2;
    [VperH, shape] = periodicizedHarmonic(grid.x, L, R0, R1);
    validation = validatePotential(VperH, shape, grid, R0);

    harmonicProblem = makeProblem(grid, H, template, 'direct harmonic');
    initial = mass * ones(N, 1) / grid.domain_length;
    fprintf('L=%g N=%d h=%.15g: solving direct harmonic ...\n', ...
        L, N, grid.h);
    harmonicResult = src.SolveGroundState1D( ...
        harmonicProblem, initial, solver);
    harmonicResult.grid = grid;
    harmonicResult.fourier = ...
        src.diagnostics.FourierTailDiagnostics(harmonicResult.rho);

    periodicProblem = makeProblem(grid, VperH, template, ...
        'harmonic with C-infinity boundary continuation');
    fprintf('L=%g N=%d: solving periodicized harmonic ...\n', L, N);
    periodicResult = src.SolveGroundState1D( ...
        periodicProblem, harmonicResult.rho, solver);
    periodicResult.grid = grid;
    periodicResult.fourier = ...
        src.diagnostics.FourierTailDiagnostics(periodicResult.rho);

    difference = periodicResult.rho - harmonicResult.rho;
    coreMask = abs(grid.x) <= R_core;
    modificationMask = abs(grid.x) >= R0;
    deltaV = VperH - H;
    harmonicEnergyPeriodicState = src.discretization.ps.Energy( ...
        periodicResult.rho, harmonicProblem);
    harmonicEnergyDirectState = src.discretization.ps.Energy( ...
        harmonicResult.rho, harmonicProblem);

    record.L = L;
    record.N = N;
    record.h = grid.h;
    record.R0 = R0;
    record.R1 = R1;
    record.interior_potential_error = validation.interior_error;
    record.endpoint_value_jump = validation.endpoint_value_jump;
    record.endpoint_derivative_jump = validation.endpoint_derivative_jump;
    record.monotonicity_pass = validation.monotonicity_pass;
    record.minimum_positive_side_difference = ...
        validation.minimum_positive_side_difference;
    record.rho_mod_max = max(periodicResult.rho(modificationMask));
    record.mass_mod = grid.h * sum(periodicResult.rho(modificationMask));
    record.p_mass_mod = grid.h * sum( ...
        p_sigma(periodicResult.rho(modificationMask)));
    record.d_global = sqrt(grid.h * sum(difference .^ 2));
    record.d_core = sqrt(grid.h * sum(difference(coreMask) .^ 2));
    record.d_core_inf = max(abs(difference(coreMask)));
    record.dE_H = abs( ...
        harmonicEnergyPeriodicState - harmonicEnergyDirectState);
    record.DeltaE_V_signed = grid.h * sum( ...
        deltaV .* p_sigma(periodicResult.rho));
    record.abs_DeltaE_V = abs(record.DeltaE_V_signed);
    record.PG_H = harmonicResult.diagnostics.final_pg_residual;
    record.KKT_H = harmonicResult.diagnostics.final_kkt_residual;
    record.exact_zeros_H = harmonicResult.diagnostics.exact_zero_count;
    record.Fourier_tail_H = harmonicResult.fourier.tail_ratio_quarter;
    record.PG_perH = periodicResult.diagnostics.final_pg_residual;
    record.KKT_perH = periodicResult.diagnostics.final_kkt_residual;
    record.exact_zeros_perH = periodicResult.diagnostics.exact_zero_count;
    record.Fourier_tail_perH = periodicResult.fourier.tail_ratio_quarter;
    record.boundary_modification_pass = ...
        record.d_core <= 0.1 * target_state_accuracy ...
        && record.dE_H <= 0.1 * target_energy_accuracy;
    record.d_core_over_target = record.d_core / target_state_accuracy;
    record.d_core_inf_over_target = ...
        record.d_core_inf / target_state_accuracy;
    record.dE_H_over_target = record.dE_H / target_energy_accuracy;
    record.DeltaE_V_over_target = ...
        record.abs_DeltaE_V / target_energy_accuracy;
    if max(record.PG_H, record.PG_perH) > 1e-8
        warning('L=%g has PG above 1e-8; do not certify its comparison.', L);
    end

    if index == 1
        records = record;
    else
        records(index) = record;
    end
    harmonicResults{index} = harmonicResult;
    periodicResults{index} = periodicResult;
    shape.validation = validation;
    potentialShapes{index} = shape;
end

hValues = [records.h];
relativeHSpread = max(abs(hValues / h_target - 1));
if relativeHSpread > 100 * eps
    error('Fixed-h comparison failed: relative h spread %.3e.', relativeHSpread);
end

% Cross-box comparisons use spectral evaluation, never real-space interp1.
xCore = linspace(-R_core, R_core, 2001).';
rhoCore = zeros(numel(xCore), caseCount);
nodalReconstructionError = zeros(caseCount, 1);
for index = 1:caseCount
    result = periodicResults{index};
    [coefficients, ~] = ...
        src.discretization.ps.FourierCoefficients(result.rho);
    rhoCore(:, index) = ...
        src.discretization.ps.EvaluateFourierSeries1D( ...
        coefficients, result.grid.L, xCore);
    nodal = src.discretization.ps.EvaluateFourierSeries1D( ...
        coefficients, result.grid.L, result.grid.x);
    nodalReconstructionError(index) = max(abs(nodal - result.rho));
end
if max(nodalReconstructionError) > 1e-12
    error('Cross-box Fourier evaluation failed nodal reconstruction.');
end
referenceCore = rhoCore(:, end);
boxComparisons = repmat(struct(), caseCount - 1, 1);
for index = 1:(caseCount - 1)
    difference = rhoCore(:, index) - referenceCore;
    boxComparisons(index).L = L_list(index);
    boxComparisons(index).reference_L = L_list(end);
    boxComparisons(index).L2_core = sqrt(trapz(xCore, difference .^ 2));
    boxComparisons(index).Linf_core = max(abs(difference));
    boxComparisons(index).pass = boxComparisons(index).L2_core ...
        <= 0.1 * target_state_accuracy;
end

recommendedL = NaN;
recommendationStatus = ...
    'no tested box is yet demonstrably large enough.';
for index = 1:caseCount
    if ~records(index).boundary_modification_pass
        continue;
    end
    if index < caseCount
        crossPass = boxComparisons(index).pass;
    else
        crossPass = true;
    end
    if crossPass
        recommendedL = records(index).L;
        if index == caseCount
            recommendationStatus = sprintf([ ...
                'L=%g is the current candidate, but a larger-box ' ...
                'confirmation may still be useful.'], recommendedL);
        else
            recommendationStatus = sprintf( ...
                'L=%g satisfies both numerical design criteria.', ...
                recommendedL);
        end
        break;
    end
end

printSensitivityTable(records);
printCrossBoxTable(boxComparisons);
printValidationSummary(records);
fprintf('\nFixed-h check: target %.15g, maximum relative spread %.3e\n', ...
    h_target, relativeHSpread);
fprintf('Maximum nodal Fourier reconstruction error: %.3e\n', ...
    max(nodalReconstructionError));
fprintf('recommended_L: %s\n', numberOrNaN(recommendedL));
fprintf('%s\n', recommendationStatus);

archive.description = ['Harmonic boundary-periodicization sensitivity ' ...
    'at fixed epsilon, sigma, and physical h.'];
archive.epsilon = epsilon;
archive.sigma = sigma;
archive.beta = beta;
archive.delta = delta;
archive.mass = mass;
archive.L_list = L_list;
archive.N_list = [records.N];
archive.h_target = h_target;
archive.R_core = R_core;
archive.target_state_accuracy = target_state_accuracy;
archive.target_energy_accuracy = target_energy_accuracy;
archive.fisher_regularization.epsilon = epsilon;
archive.fisher_regularization.s_epsilon = s_epsilon;
archive.fisher_regularization.ds_epsilon = ds_epsilon;
archive.fisher_regularization.d2s_epsilon = d2s_epsilon;
archive.fisher_regularization.label = fisher_label;
archive.potential_regularization.sigma = sigma;
archive.potential_regularization.p_sigma = p_sigma;
archive.potential_regularization.dp_sigma = dp_sigma;
archive.potential_regularization.d2p_sigma = d2p_sigma;
archive.potential_regularization.label = potential_label;
archive.records = records;
archive.box_comparisons = boxComparisons;
archive.x_core = xCore;
archive.rho_perH_core = rhoCore;
archive.nodal_reconstruction_error = nodalReconstructionError;
archive.rho_H = cellfun(@(result) result.rho, harmonicResults, ...
    'UniformOutput', false);
archive.rho_perH = cellfun(@(result) result.rho, periodicResults, ...
    'UniformOutput', false);
archive.grid = cellfun(@(result) result.grid, periodicResults, ...
    'UniformOutput', false);
archive.solver_diagnostics_H = cellfun(@(result) result.diagnostics, ...
    harmonicResults, 'UniformOutput', false);
archive.solver_diagnostics_perH = cellfun(@(result) result.diagnostics, ...
    periodicResults, 'UniformOutput', false);
archive.potential_shapes = potentialShapes;
archive.recommended_L = recommendedL;
archive.recommendation_status = recommendationStatus;
archive.solver = solver;
save(outputFile, 'archive', '-v7.3');
makeFigures(records, periodicResults, potentialShapes, ...
    xCore, rhoCore, outputFolder);

function solver = makeSolver()
config = experiments.DefaultConfig();
solver = config.solver;
solver.name = 'fista_cd';
solver.splitting = 'potential_prox';
solver.projection_name = 'semismooth';
solver.projection_tol = 1e-14;
solver.pg_tol = 1e-5;
solver.final_pg_tol = 1e-12;
solver.certification_tol = 1e-12;
solver.residual_check_interval = 10;
solver.max_iter = 200000;
solver.switch.enabled = true;
solver.switch.energy_window = 50;
solver.switch.energy_tol = 1e-12;
solver.switch.consecutive_windows = 2;
solver.switch.min_iter = 200;
solver.switch.pg_entry_tol = 1e-5;
solver.switch.max_main_iter = 20000;
solver.polish_mode = 'if_needed';
solver.polish.pg_tol = 1e-12;
solver.polish.max_iter = 20;
solver.polish.linear_solver = 'interior_pcg_schur';
solver.polish.preconditioner = 'fd_variable';
solver.polish.allow_pdas_fallback = true;
solver.display = false;
end

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end

function [V, shape] = periodicizedHarmonic(x, L, R0, R1)
shape.H = 0.5 * x .^ 2;
shape.R0 = R0;
shape.R1 = R1;
shape.Vflat = 0.5 * L ^ 2;
shape.chi = model.SmoothStepCInf( ...
    (abs(x) - R0) / (R1 - R0));
V = (1 - shape.chi) .* shape.H + shape.chi .* shape.Vflat;
shape.V = V;
end

function validation = validatePotential(V, shape, grid, R0)
interior = abs(grid.x) <= R0;
validation.interior_error = max(abs(V(interior) - shape.H(interior)));
validation.minimum_value = min(V);
positive = grid.x >= 0;
positiveDifferences = diff(V(positive));
validation.minimum_positive_side_difference = min(positiveDifferences);
monotonicityTolerance = -1e-12 * max(1, max(V));
validation.monotonicity_pass = ...
    validation.minimum_positive_side_difference >= monotonicityTolerance;
derivative = src.discretization.ps.FirstDerivative( ...
    V, src.discretization.ps.Plan1D(grid));
validation.endpoint_value_jump = abs(V(end) - V(1));
validation.endpoint_derivative_jump = abs(derivative(end) - derivative(1));
if validation.interior_error > 100 * eps(max(1, max(shape.H)))
    error('Periodic continuation changed the harmonic interior.');
end
if validation.minimum_value < -100 * eps(max(1, max(V)))
    error('Periodic continuation created a negative potential.');
end
if ~validation.monotonicity_pass
    error('Periodic continuation created an artificial potential well.');
end
end

function problem = makeProblem(grid, V, template, label)
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = V(:);
problem.beta = template.beta;
problem.delta = template.delta;
problem.mass = template.mass;
problem.fisher_regularization.epsilon = template.epsilon;
problem.fisher_regularization.s_epsilon = template.s_epsilon;
problem.fisher_regularization.ds_epsilon = template.ds_epsilon;
problem.fisher_regularization.d2s_epsilon = template.d2s_epsilon;
problem.fisher_regularization.label = template.fisher_label;
problem.potential_regularization.name = 'inline_fixed_sigma';
problem.potential_regularization.sigma = template.sigma;
problem.potential_regularization.p_sigma = template.p_sigma;
problem.potential_regularization.dp_sigma = template.dp_sigma;
problem.potential_regularization.d2p_sigma = template.d2p_sigma;
problem.potential_regularization.label = template.potential_label;
problem.potential_regularization.prox_type = 'generic_convex';
problem.potential_regularization.convexity_tol = 1e-14;
problem.potential_label = label;
problem.entropy.enabled = false;
problem.entropy.eta = 0;
end

function printSensitivityTable(records)
fprintf('\n--------------------------------------------------------------\n');
fprintf('Boundary periodicization sensitivity\n');
fprintf('--------------------------------------------------------------\n');
fprintf([' L    N       h        R0    R1   rho_mod_max mass_mod   ' ...
    'd_global   d_core     d_core_inf dE_H       DeltaE_V   PG_H/PG_perH\n']);
for index = 1:numel(records)
    r = records(index);
    fprintf(['%2g  %4d  %.7f  %4.1f  %4.1f  %.3e  %.3e  %.3e  ' ...
        '%.3e  %.3e  %.3e  %+.3e  %.2e/%.2e\n'], ...
        r.L, r.N, r.h, r.R0, r.R1, r.rho_mod_max, r.mass_mod, ...
        r.d_global, r.d_core, r.d_core_inf, r.dE_H, ...
        r.DeltaE_V_signed, r.PG_H, r.PG_perH);
    fprintf(['    target ratios: d_core %.3e, d_core_inf %.3e, ' ...
        'dE_H %.3e, |DeltaE_V| %.3e; boundary=%s\n'], ...
        r.d_core_over_target, r.d_core_inf_over_target, ...
        r.dE_H_over_target, r.DeltaE_V_over_target, ...
        passFail(r.boundary_modification_pass));
end
end

function printCrossBoxTable(comparisons)
fprintf('\n--------------------------------------------------------------\n');
fprintf('Cross-box core-state sensitivity\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' comparison       L2_core      Linf_core\n');
for index = 1:numel(comparisons)
    c = comparisons(index);
    fprintf(' L%g vs L%g     %.3e     %.3e\n', ...
        c.L, c.reference_L, c.L2_core, c.Linf_core);
end
end

function printValidationSummary(records)
fprintf('\n--------------------------------------------------------------\n');
fprintf('Potential validation summary\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' L    interior_error  value_jump   derivative_jump monotonic\n');
for index = 1:numel(records)
    r = records(index);
    fprintf('%2g    %.3e       %.3e    %.3e       %s\n', ...
        r.L, r.interior_potential_error, r.endpoint_value_jump, ...
        r.endpoint_derivative_jump, passFail(r.monotonicity_pass));
end
end

function makeFigures(records, periodicResults, shapes, xCore, rhoCore, folder)
colors = lines(numel(records));
f1 = figure('Visible', 'off', 'Name', 'Boundary potential continuation');
tiledlayout(1, numel(records));
for index = 1:numel(records)
    nexttile; hold on;
    sourceGrid = periodicResults{index}.grid;
    positive = sourceGrid.x >= 0;
    plot(sourceGrid.x(positive), shapes{index}.H(positive), '-', ...
        sourceGrid.x(positive), shapes{index}.V(positive), '--', ...
        'LineWidth', 1.1);
    xline(records(index).R0, ':', 'R0');
    xline(records(index).R1, ':', 'R1');
    grid(gca, 'on'); title(sprintf('L=%g', records(index).L));
    xlabel('x'); ylabel('V');
end
legend('harmonic', 'periodicized', 'Location', 'best');
exportgraphics(f1, fullfile(folder, '01_potential_continuation.png'), ...
    'Resolution', 180);
close(f1);

f2 = figure('Visible', 'off', 'Name', 'Core periodicized states');
tiledlayout(1, 2);
nexttile; hold on;
for index = 1:numel(records)
    plot(xCore, rhoCore(:, index), 'Color', colors(index, :), ...
        'DisplayName', sprintf('L=%g', records(index).L));
end
grid on; xlabel('x'); ylabel('rho'); title('linear scale');
legend('Location', 'best');
nexttile; hold on;
for index = 1:numel(records)
    semilogy(xCore, max(rhoCore(:, index), realmin), ...
        'Color', colors(index, :), ...
        'DisplayName', sprintf('L=%g', records(index).L));
end
grid on; xlabel('x'); ylabel('rho'); title('semilog scale');
legend('Location', 'best');
exportgraphics(f2, fullfile(folder, '02_core_states.png'), ...
    'Resolution', 180);
close(f2);

f3 = figure('Visible', 'off', 'Name', 'Boundary sensitivity');
tiledlayout(1, 2);
nexttile;
semilogy([records.L], [records.d_core], 'o-');
grid on; xlabel('L'); ylabel('core L2 state difference');
nexttile;
semilogy([records.L], [records.dE_H], 'o-', ...
    [records.L], [records.abs_DeltaE_V], 's-');
grid on; xlabel('L'); ylabel('energy diagnostic');
legend('dE_H', '|DeltaE_V|', 'Location', 'best');
exportgraphics(f3, fullfile(folder, '03_boundary_sensitivity.png'), ...
    'Resolution', 180);
close(f3);
end

function label = passFail(value)
if value
    label = 'PASS';
else
    label = 'FAIL';
end
end

function text = numberOrNaN(value)
if isfinite(value)
    text = sprintf('%g', value);
else
    text = 'NaN';
end
end
