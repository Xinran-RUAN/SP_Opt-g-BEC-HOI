%RUN_SPATIAL_OPERATOR_CONSISTENCY Diagnose termwise Fourier consistency.
%
% This run never invokes an optimizer and does not alter the production
% discretization. Every derivative, adjoint derivative, Fourier projection,
% and coefficient normalization is evaluated through the production API.
clearvars; clc;

epsilon = 1e-3;
sigma = 1e-12;
beta = 10;
delta = 10;
mass = 1;
L = 8;
N_list = [32, 64, 128, 256, 512, 1024];
N_ref_diag = 16384;
save_result = true;
save_figures = true;

s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
if any(mod(N_ref_diag, N_list) ~= 0) || N_ref_diag <= max(N_list)
    error('N_ref_diag must be a strictly finer integer multiple of every N.');
end

settings.epsilon = epsilon;
settings.sigma = sigma;
settings.beta = beta;
settings.delta = delta;
settings.mass = mass;
settings.L = L;
settings.s_epsilon = s_epsilon;
settings.ds_epsilon = ds_epsilon;
settings.d2s_epsilon = d2s_epsilon;
settings.p_sigma = p_sigma;
settings.dp_sigma = dp_sigma;
settings.d2p_sigma = d2p_sigma;

gridRef = makeGrid(L, N_ref_diag);
rhoRef = manufacturedDensity(gridRef.x, L, mass);
if min(rhoRef) <= 0
    error('The manufactured density must be strictly positive.');
end
fprintf('\nManufactured density: min %.6e, mass error %.3e, Nref %d\n', ...
    min(rhoRef), abs(gridRef.h * sum(rhoRef) - mass), N_ref_diag);

caseH = evaluateCase('harmonic', rhoRef, gridRef, N_list, settings);
caseP = evaluateCase('periodic', rhoRef, gridRef, N_list, settings);
smoothness = endpointSmoothness(settings);

fitOptions.fit_range = [8, 512];
fitOptions.relative_floor = 1e-11;
envelopes.rho_ref = src.diagnostics.FourierEnvelopeFit( ...
    rhoRef, fitOptions);
envelopes.GV_harmonic = src.diagnostics.FourierEnvelopeFit( ...
    caseH.reference_components.potential, fitOptions);
envelopes.GV_periodic = src.diagnostics.FourierEnvelopeFit( ...
    caseP.reference_components.potential, fitOptions);
envelopes.GF = src.diagnostics.FourierEnvelopeFit( ...
    caseH.reference_components.fisher, fitOptions);
envelopes.Gtotal = src.diagnostics.FourierEnvelopeFit( ...
    caseH.reference_components.total, fitOptions);

actualReference = latestActualReference(root);

printOperatorTable(caseH, 'H');
printOperatorTable(caseP, 'P');
printPotentialOversampling(caseH, 'H');
printPotentialOversampling(caseP, 'P');
printFisherOversampling(caseH);
printSmoothness(smoothness);
printEnvelopeSummary(envelopes);
printEnvelopeBands(envelopes);
if actualReference.found
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Actual p=4 finest-grid state Fourier envelope\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('file: %s\n', actualReference.file);
    printOneEnvelopeBands('actual_rho_Nref', actualReference.envelope);
    fprintf(['actual rho fit k=16:512: m=%.4f R2_alg=%.4f, ' ...
        'c=%.4e R2_exp=%.4f, count=%d\n'], ...
        actualReference.envelope.algebraic_slope, ...
        actualReference.envelope.R2_algebraic, ...
        actualReference.envelope.exponential_slope, ...
        actualReference.envelope.R2_exponential, ...
        actualReference.envelope.fit_count);
else
    warning('No p=4 Nref mesh archive was found for the actual-state fit.');
end

diagnosis = classifyDiagnosis(caseH, caseP);
fprintf('\n============================================================\n');
fprintf('AUTOMATIC SPATIAL DIAGNOSIS: %s\n%s\n', ...
    diagnosis.case, diagnosis.message);
fprintf('============================================================\n');

archive.settings = settings;
archive.N_list = N_list;
archive.N_ref_diag = N_ref_diag;
archive.manufactured_density = rhoRef;
archive.grid_ref = gridRef;
archive.harmonic = caseH;
archive.periodic = caseP;
archive.smoothness = smoothness;
archive.envelopes = envelopes;
archive.actual_reference = actualReference;
archive.diagnosis = diagnosis;

outputFolder = fullfile(root, 'results', 'diagnostics');
if (save_result || save_figures) && ~isfolder(outputFolder)
    mkdir(outputFolder);
end
if save_result
    save(fullfile(outputFolder, 'spatial_operator_consistency.mat'), ...
        'archive', '-v7.3');
end
if save_figures
    makeFigures(caseH, caseP, envelopes, actualReference, outputFolder);
end

function result = evaluateCase(caseName, rhoRef, gridRef, Nvalues, settings)
Vref = potentialValues(gridRef.x, settings.L, caseName);
problemRef = makeProblem(gridRef, Vref, settings);
referenceComponents = src.diagnostics.GradientComponents(rhoRef, problemRef);
productionGradient = src.discretization.ps.Gradient(rhoRef, problemRef);
referenceSumError = norm(productionGradient - referenceComponents.total) ...
    / max(1, norm(productionGradient));
if referenceSumError > 1e-13
    error('Reference component sum failed for %s: %.3e.', ...
        caseName, referenceSumError);
end

fieldNames = {'fisher', 'potential', 'beta', 'delta', 'total'};
numberOfN = numel(Nvalues);
errors = zeros(numberOfN, numel(fieldNames));
relativeErrors = zeros(size(errors));
potentialOversampled = zeros(numberOfN, 3);
fisherOversampled = zeros(numberOfN, 3);

for j = 1:numberOfN
    gridN = makeGrid(settings.L, Nvalues(j));
    rhoN = src.diagnostics.FourierProjectReference(rhoRef, gridRef, gridN);
    VN = potentialValues(gridN.x, settings.L, caseName);
    problemN = makeProblem(gridN, VN, settings);
    componentsN = src.diagnostics.GradientComponents(rhoN, problemN);

    for ell = 1:numel(fieldNames)
        field = fieldNames{ell};
        projectedReference = src.diagnostics.FourierProjectReference( ...
            referenceComponents.(field), gridRef, gridN);
        difference = componentsN.(field) - projectedReference;
        errors(j, ell) = discreteL2(difference, gridN.h);
        relativeErrors(j, ell) = errors(j, ell) ...
            / max(realmin, discreteL2(projectedReference, gridN.h));
    end

    for factorIndex = 1:3
        factor = [1, 2, 4];
        factor = factor(factorIndex);
        if factor == 1
            potentialOnN = componentsN.potential;
            fisherOnN = componentsN.fisher;
        else
            potentialOnN = oversampledComponent( ...
                rhoN, gridN, factor, 'potential', caseName, settings);
            fisherOnN = oversampledComponent( ...
                rhoN, gridN, factor, 'fisher', caseName, settings);
        end
        potentialReferenceN = src.diagnostics.FourierProjectReference( ...
            referenceComponents.potential, gridRef, gridN);
        fisherReferenceN = src.diagnostics.FourierProjectReference( ...
            referenceComponents.fisher, gridRef, gridN);
        potentialOversampled(j, factorIndex) = discreteL2( ...
            potentialOnN - potentialReferenceN, gridN.h);
        fisherOversampled(j, factorIndex) = discreteL2( ...
            fisherOnN - fisherReferenceN, gridN.h);
    end
end

result.name = caseName;
result.N = Nvalues(:);
result.fields = fieldNames;
result.errors = errors;
result.relative_errors = relativeErrors;
result.rates = adjacentRates(errors);
result.potential_oversampling_errors = potentialOversampled;
result.potential_oversampling_rates = adjacentRates(potentialOversampled);
result.fisher_oversampling_errors = fisherOversampled;
result.fisher_oversampling_rates = adjacentRates(fisherOversampled);
result.reference_components = referenceComponents;
result.reference_component_sum_error = referenceSumError;
end

function projected = oversampledComponent( ...
    rhoN, gridN, factor, field, caseName, settings)
gridFine = makeGrid(settings.L, factor * gridN.N);
rhoFine = src.discretization.ps.Prolong(rhoN, gridFine.N);
Vfine = potentialValues(gridFine.x, settings.L, caseName);
problemFine = makeProblem(gridFine, Vfine, settings);
componentsFine = src.diagnostics.GradientComponents(rhoFine, problemFine);
projected = src.diagnostics.FourierProjectReference( ...
    componentsFine.(field), gridFine, gridN);
end

function problem = makeProblem(grid, V, settings)
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = V(:);
problem.beta = settings.beta;
problem.delta = settings.delta;
problem.mass = settings.mass;
problem.fisher_regularization.epsilon = settings.epsilon;
problem.fisher_regularization.s_epsilon = settings.s_epsilon;
problem.fisher_regularization.ds_epsilon = settings.ds_epsilon;
problem.fisher_regularization.d2s_epsilon = settings.d2s_epsilon;
problem.fisher_regularization.label = 's_epsilon(rho) = rho + epsilon';
problem.potential_regularization.sigma = settings.sigma;
problem.potential_regularization.p_sigma = settings.p_sigma;
problem.potential_regularization.dp_sigma = settings.dp_sigma;
problem.potential_regularization.d2p_sigma = settings.d2p_sigma;
problem.potential_regularization.label = ...
    'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
problem.potential_regularization.prox_type = 'generic_convex';
end

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end

function rho = manufacturedDensity(x, L, mass)
theta = pi * x / L;
raw = 1 + 0.20 * cos(theta) + 0.10 * cos(2 * theta) ...
    + 0.05 * sin(3 * theta);
rho = (mass / (2 * L)) * raw;
end

function derivative = manufacturedDensityDerivative(x, L, mass)
theta = pi * x / L;
rawDerivative = (pi / L) * (-0.20 * sin(theta) ...
    - 0.20 * sin(2 * theta) + 0.15 * cos(3 * theta));
derivative = (mass / (2 * L)) * rawDerivative;
end

function V = potentialValues(x, L, caseName)
switch lower(caseName)
    case 'harmonic'
        V = 0.5 * x .^ 2;
    case 'periodic'
        V = (L ^ 2 / pi ^ 2) * (1 - cos(pi * x / L));
    otherwise
        error('Unknown diagnostic potential case "%s".', caseName);
end
end

function derivative = potentialDerivative(x, L, caseName)
switch lower(caseName)
    case 'harmonic'
        derivative = x;
    case 'periodic'
        derivative = (L / pi) * sin(pi * x / L);
    otherwise
        error('Unknown diagnostic potential case "%s".', caseName);
end
end

function smoothness = endpointSmoothness(settings)
x = [-settings.L; settings.L];
rho = manufacturedDensity(x, settings.L, settings.mass);
rhoX = manufacturedDensityDerivative(x, settings.L, settings.mass);
names = {'V_H', 'V_periodic'};
cases = {'harmonic', 'periodic'};
for j = 1:2
    V = potentialValues(x, settings.L, cases{j});
    Vx = potentialDerivative(x, settings.L, cases{j});
    gV = V .* settings.dp_sigma(rho);
    gVx = Vx .* settings.dp_sigma(rho) ...
        + V .* settings.d2p_sigma(rho) .* rhoX;
    smoothness.(names{j}).value_jump = V(2) - V(1);
    smoothness.(names{j}).derivative_jump = Vx(2) - Vx(1);
    smoothness.([names{j}, '_times_dp']).value_jump = gV(2) - gV(1);
    smoothness.([names{j}, '_times_dp']).derivative_jump = ...
        gVx(2) - gVx(1);
end
end

function result = latestActualReference(root)
result.found = false;
result.file = '';
result.envelope = struct();
folder = fullfile(root, 'results', 'potential_mesh');
matches = dir(fullfile(folder, ...
    'potential_mesh_inline_p_sigma_p4_Nref*.mat'));
if isempty(matches)
    return;
end
[~, newest] = max([matches.datenum]);
file = fullfile(matches(newest).folder, matches(newest).name);
loaded = load(file, 'rho');
if ~isfield(loaded, 'rho') || isempty(loaded.rho)
    return;
end
rho = loaded.rho;
if iscell(rho)
    rho = rho{end};
end
fitOptions.fit_range = [16, 512];
fitOptions.relative_floor = 100 * eps;
result.found = true;
result.file = file;
result.envelope = src.diagnostics.FourierEnvelopeFit(rho, fitOptions);
end

function normValue = discreteL2(values, h)
normValue = sqrt(h * sum(abs(values) .^ 2));
end

function rates = adjacentRates(errors)
rates = NaN(size(errors));
for j = 2:size(errors, 1)
    valid = errors(j - 1, :) > 0 & errors(j, :) > 0;
    rates(j, valid) = log(errors(j - 1, valid) ...
        ./ errors(j, valid)) / log(2);
end
end

function printOperatorTable(result, shortName)
fprintf('\n------------------------------------------------------------\n');
fprintf('Manufactured operator consistency: Case %s (%s)\n', ...
    shortName, result.name);
fprintf('------------------------------------------------------------\n');
fprintf(' N        eF           eV          ebeta        edelta       etotal\n');
for j = 1:numel(result.N)
    fprintf('%4d  %.3e  %.3e  %.3e  %.3e  %.3e\n', ...
        result.N(j), result.errors(j, :));
end
fprintf('relative errors\n');
for j = 1:numel(result.N)
    fprintf('%4d  %.3e  %.3e  %.3e  %.3e  %.3e\n', ...
        result.N(j), result.relative_errors(j, :));
end
fprintf('adjacent rates log2(e_N/e_2N)\n');
for j = 2:numel(result.N)
    fprintf('%4d  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n', ...
        result.N(j), result.rates(j, :));
end
end

function printPotentialOversampling(result, shortName)
fprintf('\n------------------------------------------------------------\n');
fprintf('Potential product: Case %s\n', shortName);
fprintf('------------------------------------------------------------\n');
fprintf(' N        native         2x             4x\n');
for j = 1:numel(result.N)
    fprintf('%4d  %.3e  %.3e  %.3e\n', result.N(j), ...
        result.potential_oversampling_errors(j, :));
end
fprintf('rates\n');
for j = 2:numel(result.N)
    fprintf('%4d  %8.3f  %8.3f  %8.3f\n', result.N(j), ...
        result.potential_oversampling_rates(j, :));
end
end

function printFisherOversampling(result)
fprintf('\n------------------------------------------------------------\n');
fprintf('Fisher product oversampling (same for H/P)\n');
fprintf('------------------------------------------------------------\n');
fprintf(' N        native         2x             4x\n');
for j = 1:numel(result.N)
    fprintf('%4d  %.3e  %.3e  %.3e\n', result.N(j), ...
        result.fisher_oversampling_errors(j, :));
end
end

function printSmoothness(smoothness)
fprintf('\n------------------------------------------------------------\n');
fprintf('Periodic smoothness\n');
fprintf('------------------------------------------------------------\n');
fprintf(' quantity                 value_jump       derivative_jump\n');
names = {'V_H', 'V_periodic', 'V_H_times_dp', 'V_periodic_times_dp'};
labels = {'V_H', 'V_per', 'V_H * p''(rho)', 'V_per * p''(rho)'};
% labels above denote p-prime in plain console text; avoid quote parsing.
labels{3} = 'V_H * dp_sigma(rho)';
labels{4} = 'V_per * dp_sigma(rho)';
for j = 1:numel(names)
    entry = smoothness.(names{j});
    fprintf(' %-23s  %+13.6e  %+16.6e\n', labels{j}, ...
        entry.value_jump, entry.derivative_jump);
end
end

function printEnvelopeSummary(envelopes)
fprintf('\n------------------------------------------------------------\n');
fprintf('Fourier envelope fits (descriptive only)\n');
fprintf('------------------------------------------------------------\n');
fprintf(' quantity          alg_slope   R2_alg    exp_slope    R2_exp  count\n');
names = fieldnames(envelopes);
for j = 1:numel(names)
    entry = envelopes.(names{j});
    fprintf(' %-15s  %9.4f  %8.4f  %11.4e  %8.4f  %5d\n', ...
        names{j}, entry.algebraic_slope, entry.R2_algebraic, ...
        entry.exponential_slope, entry.R2_exponential, entry.fit_count);
end
end

function printEnvelopeBands(envelopes)
names = fieldnames(envelopes);
for j = 1:numel(names)
    printOneEnvelopeBands(names{j}, envelopes.(names{j}));
end
end

function printOneEnvelopeBands(label, envelope)
fprintf('\nFourier bands: %s\n', label);
fprintf(' |k| band       max |hat f_k|      rms |hat f_k|\n');
for j = 1:size(envelope.bands, 1)
    fprintf(' %3d:%-3d       %.6e        %.6e\n', ...
        envelope.bands(j, 1), envelope.bands(j, 2), ...
        envelope.band_max(j), envelope.band_rms(j));
end
end

function diagnosis = classifyDiagnosis(caseH, caseP)
betaIndex = find(strcmp(caseH.fields, 'beta'));
deltaIndex = find(strcmp(caseH.fields, 'delta'));
fisherIndex = find(strcmp(caseH.fields, 'fisher'));
potentialIndex = find(strcmp(caseH.fields, 'potential'));
linearSanity = max(caseH.errors(:, [betaIndex, deltaIndex]), [], 'all') ...
    <= 1e-10;
fisherRapid = caseH.errors(end, fisherIndex) <= 1e-10;
harmonicAlgebraic = caseH.errors(end, potentialIndex) > 1e-10;
periodicRapid = caseP.errors(end, potentialIndex) <= 1e-10;
nativeFinal = caseH.potential_oversampling_errors(end, 1);
fourXFinal = caseH.potential_oversampling_errors(end, 3);
oversamplingRestores = fourXFinal <= 1e-10 ...
    || fourXFinal <= 1e-4 * nativeFinal;

if ~linearSanity
    diagnosis.case = 'CASE E';
    diagnosis.message = ['Fourier indexing/projection/normalization is ' ...
        'inconsistent; stop higher-level interpretation.'];
elseif harmonicAlgebraic && oversamplingRestores
    diagnosis.case = 'CASE B';
    diagnosis.message = ['Native pseudospectral product aliasing is the ' ...
        'dominant obstruction.'];
elseif ~fisherRapid && caseH.fisher_oversampling_errors(end, 3) > 1e-10
    diagnosis.case = 'CASE C';
    diagnosis.message = ['The Fisher term remains algebraic after ' ...
        'oversampling; inspect Fisher discretization/regularization.'];
elseif fisherRapid && harmonicAlgebraic && periodicRapid
    diagnosis.case = 'CASE A';
    diagnosis.message = ['The non-smooth periodic extension of harmonic ' ...
        'V is the dominant spectral-consistency obstruction.'];
else
    diagnosis.case = 'CASE D';
    diagnosis.message = ['The manufactured operator is spectrally ' ...
        'consistent; investigate minimizer regularity, boundary ' ...
        'compatibility, or reference comparison.'];
end
diagnosis.linear_sanity = linearSanity;
diagnosis.fisher_rapid = fisherRapid;
diagnosis.harmonic_potential_algebraic = harmonicAlgebraic;
diagnosis.periodic_potential_rapid = periodicRapid;
diagnosis.oversampling_restores = oversamplingRestores;
end

function makeFigures(caseH, caseP, envelopes, actualReference, folder)
N = caseH.N;
labels = {'Fisher', 'potential', 'beta', 'delta', 'total'};
f1 = figure('Visible', 'off', 'Name', 'Operator consistency');
tiledlayout(1, 2);
nexttile;
loglog(N, max(caseH.errors, realmin), 'o-', 'LineWidth', 1.1);
grid on; xlabel('N'); ylabel('L2 operator error');
title('harmonic V'); legend(labels, 'Location', 'best');
nexttile;
loglog(N, max(caseP.errors, realmin), 'o-', 'LineWidth', 1.1);
grid on; xlabel('N'); ylabel('L2 operator error');
title('smooth periodic V'); legend(labels, 'Location', 'best');
exportgraphics(f1, fullfile(folder, 'spatial_operator_components.png'), ...
    'Resolution', 180);
close(f1);

f2 = figure('Visible', 'off', 'Name', 'Potential oversampling');
tiledlayout(1, 2);
nexttile;
loglog(N, caseH.potential_oversampling_errors, 'o-', 'LineWidth', 1.1);
grid on; xlabel('N'); ylabel('L2 potential error'); title('harmonic V');
legend('native', '2x', '4x', 'Location', 'best');
nexttile;
loglog(N, max(caseP.potential_oversampling_errors, realmin), ...
    'o-', 'LineWidth', 1.1);
grid on; xlabel('N'); ylabel('L2 potential error'); title('periodic V');
legend('native', '2x', '4x', 'Location', 'best');
exportgraphics(f2, fullfile(folder, 'spatial_operator_oversampling.png'), ...
    'Resolution', 180);
close(f2);

f3 = figure('Visible', 'off', 'Name', 'Reference Fourier envelopes');
hold on;
names = fieldnames(envelopes);
for j = 1:numel(names)
    entry = envelopes.(names{j});
    positive = entry.modes > 0;
    semilogy(entry.modes(positive), entry.amplitudes(positive), ...
        'DisplayName', names{j});
end
if actualReference.found
    entry = actualReference.envelope;
    positive = entry.modes > 0;
    semilogy(entry.modes(positive), entry.amplitudes(positive), ...
        'DisplayName', 'actual rho Nref');
end
xlim([1, 512]); grid on; xlabel('|k|'); ylabel('|hat f_k|');
legend('Location', 'best');
exportgraphics(f3, fullfile(folder, 'spatial_operator_envelopes.png'), ...
    'Resolution', 180);
close(f3);
end
