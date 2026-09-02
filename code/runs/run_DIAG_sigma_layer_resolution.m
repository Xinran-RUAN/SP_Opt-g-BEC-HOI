%RUN_DIAG_SIGMA_LAYER_RESOLUTION Small-sigma transition-layer diagnosis.
%
% This diagnostic constructs temporary problems but never changes the
% production potential, solver, regularization, or mesh-refinement code.
clearvars; clc;

epsilon = 1e-3;
sigma_list = [1e-3, 1e-6, 1e-9, 1e-12];
beta = 10;
delta = 10;
mass = 1;
L = 8;
N_fixed = 2048;
mesh_N = [128, 256, 512, 1024, 2048];
N_ref = 4096;
periodic_mesh_sigma_list = [1e-6, 1e-9, 1e-12];
save_figures = true;
save_result = true;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
solver = diagnosticSolver();
checkpointFile = fullfile(root, 'results', 'diagnostics', ...
    'sigma_layer_fixed_checkpoint.mat');

fprintf('\nSmall-sigma transition-layer diagnostic\n');
fprintf('epsilon %.3e, L %g, fixed N %d, mesh reference N %d\n', ...
    epsilon, L, N_fixed, N_ref);

% The sigma=1e-12 harmonic N=2048 state already exists in the production mesh
% archive. It is read rather than recomputed.
fixed.H = cell(numel(sigma_list), 1);
fixed.P = cell(numel(sigma_list), 1);
fixed.H{4} = loadHarmonicSmallSigma(root, N_fixed, sigma_list(4), ...
    epsilon, beta, delta, mass, L, solver);

% Linear baseline and the first three sigma cases are checkpointed states.
if isfile(checkpointFile)
    checkpointData = load(checkpointFile, 'checkpoint');
    linearProblem = makeProblem(N_fixed, 'harmonic', 0, epsilon, ...
        beta, delta, mass, L);
    linearResult = compactCachedResult( ...
        checkpointData.checkpoint.linear_rho, linearProblem, solver);
    for index = 1:3
        problem = makeProblem(N_fixed, 'harmonic', sigma_list(index), ...
            epsilon, beta, delta, mass, L);
        fixed.H{index} = compactCachedResult( ...
            checkpointData.checkpoint.harmonic_rho{index}, problem, solver);
    end
    fprintf('Reused N=2048 linear and harmonic sigma checkpoint.\n');
else
    linearProblem = makeProblem(N_fixed, 'harmonic', 0, epsilon, ...
        beta, delta, mass, L);
    linearInitial = fixed.H{4}.rho;
    fprintf('Solving diagnostic linear baseline at N=%d...\n', N_fixed);
    linearResult = src.SolveGroundState1D( ...
        linearProblem, linearInitial, solver);
    linearResult.grid = linearProblem.grid;
    for index = 1:3
        sigma = sigma_list(index);
        problem = makeProblem(N_fixed, 'harmonic', sigma, epsilon, ...
            beta, delta, mass, L);
        initial = loadHarmonicSigmaWarm(root, index, problem.grid, mass);
        fprintf('Solving H, sigma=%.0e, N=%d...\n', sigma, N_fixed);
        fixed.H{index} = src.SolveGroundState1D(problem, initial, solver);
        fixed.H{index}.grid = problem.grid;
    end
    checkpoint.linear_rho = linearResult.rho;
    checkpoint.harmonic_rho = cell(3, 1);
    for index = 1:3
        checkpoint.harmonic_rho{index} = fixed.H{index}.rho;
    end
    checkpoint.harmonic_rho{4} = fixed.H{4}.rho;
    checkpointFolder = fileparts(checkpointFile);
    if ~isfolder(checkpointFolder)
        mkdir(checkpointFolder);
    end
    save(checkpointFile, 'checkpoint', '-v7.3');
end
freeBoundary = estimateFreeBoundary(linearResult.rho, linearProblem.grid);
fprintf('linear free-boundary estimate: x_FB=%.12g (%s)\n', ...
    freeBoundary.x_right, freeBoundary.status);

% Periodic potential: continuation from smaller to larger sigma.
checkpointData = load(checkpointFile, 'checkpoint');
if isfield(checkpointData.checkpoint, 'periodic_rho')
    for index = 2:4
        problem = makeProblem(N_fixed, 'periodic', sigma_list(index), ...
            epsilon, beta, delta, mass, L);
        fixed.P{index} = compactCachedResult( ...
            checkpointData.checkpoint.periodic_rho{index}, problem, solver);
    end
    fprintf('Reused N=2048 periodic sigma checkpoint.\n');
else
    periodicOrder = [4, 3, 2];
    previousPeriodic = [];
    for orderIndex = 1:numel(periodicOrder)
        index = periodicOrder(orderIndex);
        sigma = sigma_list(index);
        problem = makeProblem(N_fixed, 'periodic', sigma, epsilon, ...
            beta, delta, mass, L);
        if index == 4
            initial = loadPeriodicP4Warm(root, problem.grid, mass);
        else
            initial = previousPeriodic.rho;
        end
        fprintf('Solving P, sigma=%.0e, N=%d...\n', sigma, N_fixed);
        fixed.P{index} = src.SolveGroundState1D(problem, initial, solver);
        fixed.P{index}.grid = problem.grid;
        previousPeriodic = fixed.P{index};
    end
    checkpoint = checkpointData.checkpoint;
    checkpoint.periodic_rho = cell(4, 1);
    for index = 2:4
        checkpoint.periodic_rho{index} = fixed.P{index}.rho;
    end
    save(checkpointFile, 'checkpoint', '-v7.3');
end

for index = 1:4
    fixed.analysisH{index, 1} = analyzeResult( ...
        fixed.H{index}, 'harmonic', sigma_list(index), ...
        freeBoundary.x_right);
end
for index = 2:4
    fixed.analysisP{index, 1} = analyzeResult( ...
        fixed.P{index}, 'periodic', sigma_list(index), ...
        freeBoundary.x_right);
end

scaling.primaryH = fitWidthScaling( ...
    sigma_list, cellfun(@(a) a.layer.primary.width_mean, fixed.analysisH));
scaling.secondaryH = fitWidthScaling( ...
    sigma_list, cellfun(@(a) a.layer.secondary.width_mean, fixed.analysisH));
scaling.slopeH = fitWidthScaling( ...
    sigma_list(2:4), cellfun(@(a) a.layer.width_slope, fixed.analysisH(2:4)));
scaling.primaryP = fitWidthScaling( ...
    sigma_list(2:4), cellfun(@(a) a.layer.primary.width_mean, fixed.analysisP(2:4)));
scaling.secondaryP = fitWidthScaling( ...
    sigma_list(2:4), cellfun(@(a) a.layer.secondary.width_mean, fixed.analysisP(2:4)));
scaling.slopeP = fitWidthScaling( ...
    sigma_list(2:4), cellfun(@(a) a.layer.width_slope, fixed.analysisP(2:4)));

printFixedTable(fixed.analysisH, fixed.analysisP, freeBoundary);
printScaling(scaling);
printRhoBands(fixed.analysisH, fixed.analysisP);
printGradientEnvelopeSummary(fixed.analysisH, fixed.analysisP);

% Second stage: fixed sigma on smooth periodic V. Existing sigma=1e-12 states at
% N<=1024 are reused from the previous diagnostic archive.
mesh = cell(numel(periodic_mesh_sigma_list), 1);
for sigmaIndex = 1:numel(periodic_mesh_sigma_list)
    sigma = periodic_mesh_sigma_list(sigmaIndex);
    fixedIndex = find(sigma_list == sigma, 1);
    fprintf('\nPeriodic mesh sweep sigma=%.0e\n', sigma);
    states = cell(numel(mesh_N) + 1, 1);
    previous = [];
    for nIndex = 1:numel(mesh_N)
        N = mesh_N(nIndex);
        problem = makeProblem(N, 'periodic', sigma, epsilon, ...
            beta, delta, mass, L);
        cached = [];
        if sigma == 1e-12 && N <= 1024
            cached = loadPeriodicP4State(root, N, problem, solver);
        elseif N == N_fixed
            cached = fixed.P{fixedIndex};
        end
        if ~isempty(cached)
            states{nIndex} = cached;
            fprintf('  N=%d reused existing state, PG %.3e\n', ...
                N, cached.diagnostics.final_pg_residual);
        else
            if isempty(previous)
                source = fixed.P{fixedIndex};
                initial = src.diagnostics.FourierProjectReference( ...
                    source.rho, source.grid, problem.grid);
            else
                initial = src.discretization.ps.Prolong(previous.rho, N);
            end
            initial = src.constraints.ProjectPositiveConservative( ...
                initial, mass, problem.grid.h, 1e-14);
            states{nIndex} = src.SolveGroundState1D(problem, initial, solver);
            states{nIndex}.grid = problem.grid;
            fprintf('  N=%d solved, PG %.3e\n', ...
                N, states{nIndex}.diagnostics.final_pg_residual);
        end
        previous = states{nIndex};
    end
    problemRef = makeProblem(N_ref, 'periodic', sigma, epsilon, ...
        beta, delta, mass, L);
    initialRef = src.discretization.ps.Prolong(previous.rho, N_ref);
    initialRef = src.constraints.ProjectPositiveConservative( ...
        initialRef, mass, problemRef.grid.h, 1e-14);
    fprintf('  solving comparison state Nref=%d...\n', N_ref);
    states{end} = src.SolveGroundState1D(problemRef, initialRef, solver);
    states{end}.grid = problemRef.grid;
    mesh{sigmaIndex} = makeMeshRecords(states, mesh_N, sigma);
    printMeshTable(mesh{sigmaIndex});
end

empiricalSigmaMinimum = inferSigmaMinimum( ...
    scaling.slopeP, [32, 64, 128, 256, 512, 1024, 2048], L, 8);
printSigmaMinimum(empiricalSigmaMinimum);
diagnosis = classifyLayerDiagnosis(fixed, mesh, scaling);
fprintf('\n============================================================\n');
fprintf('SIGMA-LAYER DIAGNOSIS: %s\n%s\n', ...
    diagnosis.case, diagnosis.message);
fprintf('============================================================\n');

archive.settings.epsilon = epsilon;
archive.settings.sigma_list = sigma_list;
archive.settings.beta = beta;
archive.settings.delta = delta;
archive.settings.mass = mass;
archive.settings.L = L;
archive.settings.N_fixed = N_fixed;
archive.settings.mesh_N = mesh_N;
archive.settings.N_ref = N_ref;
archive.linear_free_boundary = freeBoundary;
archive.fixed.analysisH = fixed.analysisH;
archive.fixed.analysisP = fixed.analysisP;
archive.fixed.rhoH = cellfun(@(s) s.rho, fixed.H, 'UniformOutput', false);
archive.fixed.rhoP = cell(size(fixed.P));
for index = 2:4
    archive.fixed.rhoP{index} = fixed.P{index}.rho;
end
archive.mesh = mesh;
archive.scaling = scaling;
archive.empirical_sigma_minimum = empiricalSigmaMinimum;
archive.diagnosis = diagnosis;

outputFolder = fullfile(root, 'results', 'diagnostics');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
if save_result
    save(fullfile(outputFolder, 'sigma_layer_resolution.mat'), ...
        'archive', '-v7.3');
end
if save_figures
    makeFigures(fixed, mesh, outputFolder, sigma_list);
end

function solver = diagnosticSolver()
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
solver.potential_prox.mass_tol = 1e-14;
solver.potential_prox.inner_tol = 1e-14;
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

function problem = makeProblem(N, caseName, sigma, epsilon, ...
    beta, delta, mass, L)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
switch lower(caseName)
    case 'harmonic'
        problem.V = 0.5 * grid.x .^ 2;
    case 'periodic'
        problem.V = (L ^ 2 / pi ^ 2) * (1 - cos(pi * grid.x / L));
    otherwise
        error('Unknown potential case "%s".', caseName);
end
problem.beta = beta;
problem.delta = delta;
problem.mass = mass;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(rho) rho + epsilon;
problem.fisher_regularization.ds_epsilon = @(rho) ones(size(rho));
problem.fisher_regularization.d2s_epsilon = @(rho) zeros(size(rho));
problem.fisher_regularization.label = 's_epsilon(rho) = rho + epsilon';
if sigma == 0
    problem.potential_regularization = src.potential.MakeLinear();
else
    problem.potential_regularization.sigma = sigma;
    problem.potential_regularization.name = 'inline_p_sigma';
    problem.potential_regularization.p_sigma = @(rho) ...
        rho .^ 2 ./ (hypot(rho, sigma) + sigma);
    problem.potential_regularization.dp_sigma = @(rho) ...
        rho ./ hypot(rho, sigma);
    problem.potential_regularization.d2p_sigma = @(rho) ...
        (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
    problem.potential_regularization.label = ...
        'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
    problem.potential_regularization.prox_type = 'generic_convex';
end
problem.entropy.enabled = false;
problem.entropy.eta = 0;
end

function result = compactCachedResult(rho, problem, solver)
rho = rho(:);
gradient = src.discretization.ps.Gradient(rho, problem);
[fullPg, ~] = src.solvers.FullGradientMapping( ...
    rho, gradient, problem, solver);
result.rho = rho;
result.grid = problem.grid;
result.problem = problem;
result.target_energy = src.discretization.ps.Energy(rho, problem);
result.baseline_energy = src.discretization.ps.BaselineEnergy(rho, problem);
result.diagnostics.final_pg_residual = fullPg;
result.diagnostics.main_iterations = NaN;
result.diagnostics.polish_iterations = NaN;
result.diagnostics.source = 'existing MAT state';
result.solver = solver;
end

function result = loadHarmonicSmallSigma(root, N, sigma, epsilon, ...
    beta, delta, mass, L, solver)
file = fullfile(root, 'results', 'potential_mesh', ...
    'potential_mesh_inline_p_sigma_p4_Nref2048_eps_1em03.mat');
data = load(file, 'rho', 'grid');
Nvalues = cellfun(@(g) g.N, data.grid);
index = find(Nvalues == N, 1);
if isempty(index)
    error('Existing harmonic p=4 state at N=%d was not found.', N);
end
problem = makeProblem(N, 'harmonic', sigma, epsilon, beta, delta, mass, L);
result = compactCachedResult(data.rho{index}, problem, solver);
end

function initial = loadHarmonicSigmaWarm(root, legacyIndex, gridFine, mass)
file = fullfile(root, 'results', 'potential_sigma', ...
    'potential_sigma_N512_eps_1em03.mat');
data = load(file, 'rho');
initial = src.discretization.ps.Prolong(data.rho{legacyIndex}, gridFine.N);
initial = src.constraints.ProjectPositiveConservative( ...
    initial, mass, gridFine.h, 1e-14);
end

function initial = loadPeriodicP4Warm(root, gridFine, mass)
file = fullfile(root, 'results', 'diagnostics', ...
    'DIAG_periodic_potential_mesh.mat');
data = load(file, 'archive');
initial = data.archive.rho{end};
initial = src.discretization.ps.Prolong(initial, gridFine.N);
initial = src.constraints.ProjectPositiveConservative( ...
    initial, mass, gridFine.h, 1e-14);
end

function result = loadPeriodicP4State(root, N, problem, solver)
file = fullfile(root, 'results', 'diagnostics', ...
    'DIAG_periodic_potential_mesh.mat');
if ~isfile(file)
    result = [];
    return;
end
data = load(file, 'archive');
Nvalues = cellfun(@(g) g.N, data.archive.grid);
index = find(Nvalues == N, 1);
if isempty(index)
    result = [];
else
    result = compactCachedResult(data.archive.rho{index}, problem, solver);
end
end

function freeBoundary = estimateFreeBoundary(rho, grid)
rho = rho(:);
center = find(grid.x >= 0, 1, 'first');
right = center:grid.N;
firstZeroRelative = find(rho(right) == 0, 1, 'first');
freeBoundary.status = 'exact zero transition found';
if isempty(firstZeroRelative) || firstZeroRelative == 1
    freeBoundary.x_right = NaN;
    freeBoundary.status = 'no positive-to-exact-zero transition';
    return;
end
firstZero = right(firstZeroRelative);
lastPositive = firstZero - 1;
fitIndices = max(center, lastPositive - 3):lastPositive;
degree = min(2, numel(fitIndices) - 1);
polynomial = polyfit(grid.x(fitIndices), rho(fitIndices), degree);
candidateRoots = roots(polynomial);
candidateRoots = real(candidateRoots(abs(imag(candidateRoots)) ...
    <= 1e-10 * max(1, abs(real(candidateRoots)))));
lower = grid.x(lastPositive);
upper = grid.x(firstZero) + 2 * grid.h;
candidateRoots = candidateRoots(candidateRoots >= lower ...
    & candidateRoots <= upper);
if isempty(candidateRoots)
    x1 = grid.x(lastPositive - 1);
    x2 = grid.x(lastPositive);
    y1 = rho(lastPositive - 1);
    y2 = rho(lastPositive);
    freeBoundary.x_right = x2 - y2 * (x2 - x1) / (y2 - y1);
    freeBoundary.status = 'linear extrapolation fallback';
else
    [~, closest] = min(abs(candidateRoots - grid.x(firstZero)));
    freeBoundary.x_right = candidateRoots(closest);
    freeBoundary.status = 'quadratic local extrapolation';
end
freeBoundary.last_positive_x = grid.x(lastPositive);
freeBoundary.first_zero_x = grid.x(firstZero);
freeBoundary.fit_indices = fitIndices;
end

function analysis = analyzeResult(result, caseName, sigma, xFreeBoundary)
problem = result.problem;
potential = problem.potential_regularization;
analysis.case_name = caseName;
analysis.sigma = sigma;
analysis.N = problem.grid.N;
analysis.h = problem.grid.h;
analysis.layer = src.diagnostics.TransitionLayerDiagnostics( ...
    result.rho, problem.grid, problem.plan, potential);
endpoint = endpointData(caseName, problem.grid.L);
analysis.endpoint = src.diagnostics.ActualEndpointCompatibility( ...
    result.rho, problem, endpoint);
analysis.far_field = src.diagnostics.FarFieldDiagnostics( ...
    result.rho, problem.grid, problem.plan);
analysis.components = src.diagnostics.GradientComponents(result.rho, problem);
fitOptions.fit_range = [8, min(512, problem.grid.N / 2 - 1)];
fitOptions.relative_floor = 100 * eps;
fitOptions.bands = [8, 16; 16, 32; 32, 64; 64, 128; ...
    128, 256; 256, 512; 512, 1024];
analysis.envelopes.rho = src.diagnostics.FourierEnvelopeFit( ...
    result.rho, fitOptions);
names = {'fisher', 'potential', 'beta', 'delta', 'total'};
for j = 1:numel(names)
    analysis.envelopes.(names{j}) = ...
        src.diagnostics.FourierEnvelopeFit( ...
        analysis.components.(names{j}), fitOptions);
end
analysis.min_density = min(result.rho);
analysis.final_pg = result.diagnostics.final_pg_residual;
analysis.target_energy = result.target_energy;
analysis.baseline_energy = result.baseline_energy;
analysis.x_free_boundary = xFreeBoundary;
analysis.x_layer_minus_free_boundary = ...
    analysis.layer.primary.center_right - xFreeBoundary;
analysis.x_q09_minus_free_boundary = ...
    analysis.layer.primary.x_high_right - xFreeBoundary;
end

function data = endpointData(caseName, L)
switch lower(caseName)
    case 'harmonic'
        data.V = [0.5 * L ^ 2, 0.5 * L ^ 2];
        data.dV = [-L, L];
    case 'periodic'
        boundaryValue = 2 * L ^ 2 / pi ^ 2;
        data.V = [boundaryValue, boundaryValue];
        data.dV = [0, 0];
end
end

function fit = fitWidthScaling(sigma, width)
sigma = sigma(:);
width = width(:);
valid = isfinite(sigma) & sigma > 0 & isfinite(width) & width > 0;
fit.count = nnz(valid);
fit.alpha = NaN;
fit.C = NaN;
fit.R2 = NaN;
if fit.count < 2
    return;
end
coefficients = polyfit(log(sigma(valid)), log(width(valid)), 1);
prediction = polyval(coefficients, log(sigma(valid)));
observed = log(width(valid));
fit.alpha = coefficients(1);
fit.C = exp(coefficients(2));
fit.R2 = 1 - sum((observed - prediction) .^ 2) ...
    / sum((observed - mean(observed)) .^ 2);
fit.sigma = sigma(valid);
fit.width = width(valid);
end

function records = makeMeshRecords(states, Nvalues, sigma)
reference = states{end};
records = repmat(struct(), numel(Nvalues), 1);
for j = 1:numel(Nvalues)
    result = states{j};
    comparison = src.diagnostics.SpectralStateComparison( ...
        result.rho, result.grid, reference.rho, reference.grid);
    analysis = analyzeResult(result, 'periodic', sigma, NaN);
    records(j).sigma = sigma;
    records(j).N = Nvalues(j);
    records(j).h = result.grid.h;
    records(j).width = analysis.layer.primary.width_mean;
    records(j).points = analysis.layer.primary.points_per_layer;
    records(j).censored_points_lower_bound = ...
        analysis.layer.primary.censored_points_lower_bound;
    records(j).width_slope = analysis.layer.width_slope;
    records(j).slope_points = analysis.layer.points_per_slope_width;
    records(j).secondary_width = analysis.layer.secondary.width_mean;
    records(j).secondary_points = analysis.layer.secondary.points_per_layer;
    records(j).resolved_L2 = comparison.resolved_L2_error;
    records(j).tail_L2 = comparison.reference_tail_L2;
    records(j).total_L2 = comparison.total_spectral_L2_error;
    records(j).m_rho = analysis.envelopes.rho.algebraic_slope;
    records(j).R2_alg = analysis.envelopes.rho.R2_algebraic;
    records(j).c_rho = analysis.envelopes.rho.exponential_slope;
    records(j).R2_exp = analysis.envelopes.rho.R2_exponential;
    records(j).final_pg = result.diagnostics.final_pg_residual;
end
for j = 1:numel(records)
    if j == 1
        records(j).state_rate = NaN;
    else
        records(j).state_rate = log(records(j - 1).total_L2 ...
            / records(j).total_L2) / log(2);
    end
end
rates = [records.state_rate];
onset = find(rates >= 3.5, 1, 'first');
if isempty(onset)
    onset = NaN;
else
    onset = Nvalues(onset);
end
for j = 1:numel(records)
    records(j).rapid_decay_onset_N = onset;
end
end

function printFixedTable(analysisH, analysisP, freeBoundary)
fprintf('\n------------------------------------------------------------\n');
fprintf('Transition layer at N=2048\n');
fprintf('------------------------------------------------------------\n');
fprintf(['pot sigma    width01-09 pts01-09 lowerPts slopeW slopePts ' ...
    'nodes minrho tailmass m_rho R2a c_rho R2e dGVjump PG\n']);
for index = 1:4
    printFixedRow(analysisH{index}, 'H');
end
for index = 2:4
    printFixedRow(analysisP{index}, 'P');
end
fprintf('linear diagnostic x_FB = %.12g (%s)\n', ...
    freeBoundary.x_right, freeBoundary.status);
fprintf(['Note: NaN primary widths mean q=0.1 was not reached inside ' ...
    '[-L,L]; lowerPts is a censored lower bound, not a width estimate.\n']);
end

function printFixedRow(a, potentialLabel)
fprintf([' %s  %.0e  %.3e %8.2f %8.2f %.3e %8.3g %5d ' ...
    '%.2e %.2e %6.2f %.3f %.2e %.3f %+.2e %.2e\n'], ...
    potentialLabel, a.sigma, a.layer.primary.width_mean, ...
    a.layer.primary.points_per_layer, ...
    a.layer.primary.censored_points_lower_bound, ...
    a.layer.width_slope, a.layer.points_per_slope_width, ...
    a.layer.primary.n_nodes_layer, a.min_density, ...
    a.far_field.tail_mass, a.envelopes.rho.algebraic_slope, ...
    a.envelopes.rho.R2_algebraic, ...
    a.envelopes.rho.exponential_slope, ...
    a.envelopes.rho.R2_exponential, ...
    a.endpoint.GV_derivative_jump, a.final_pg);
fprintf('      q09-xFB=%+.3e, primary=%s, secondary width=%.3e pts=%.2f\n', ...
    a.x_q09_minus_free_boundary, a.layer.primary.reason, ...
    a.layer.secondary.width_mean, a.layer.secondary.points_per_layer);
end

function printScaling(scaling)
fprintf('\n------------------------------------------------------------\n');
fprintf('Width scaling width = C sigma^alpha\n');
fprintf('------------------------------------------------------------\n');
names = fieldnames(scaling);
for j = 1:numel(names)
    fit = scaling.(names{j});
    fprintf('%-12s alpha=%+.6f R2=%.6f C=%.3e count=%d\n', ...
        names{j}, fit.alpha, fit.R2, fit.C, fit.count);
end
end

function printRhoBands(analysisH, analysisP)
fprintf('\n------------------------------------------------------------\n');
fprintf('Actual rho Fourier bands at N=2048\n');
fprintf('------------------------------------------------------------\n');
for index = 1:4
    printOneRhoBands(analysisH{index}, 'H');
end
for index = 2:4
    printOneRhoBands(analysisP{index}, 'P');
end
end

function printOneRhoBands(a, label)
e = a.envelopes.rho;
fprintf('%s sigma=%.0e: m=%.4f R2a=%.4f c=%.3e R2e=%.4f\n', ...
    label, a.sigma, e.algebraic_slope, e.R2_algebraic, ...
    e.exponential_slope, e.R2_exponential);
for j = 1:size(e.bands, 1)
    fprintf('  %4d:%-4d max %.3e rms %.3e\n', ...
        e.bands(j, 1), e.bands(j, 2), e.band_max(j), e.band_rms(j));
end
end

function printGradientEnvelopeSummary(analysisH, analysisP)
fprintf('\n------------------------------------------------------------\n');
fprintf('Euler-gradient component envelope fits at N=2048\n');
fprintf('------------------------------------------------------------\n');
fprintf('pot sigma term       m_alg R2_alg c_exp R2_exp\n');
for index = 1:4
    printGradientRows(analysisH{index}, 'H');
end
for index = 2:4
    printGradientRows(analysisP{index}, 'P');
end
end

function printGradientRows(a, label)
names = {'fisher', 'potential', 'beta', 'delta', 'total'};
for j = 1:numel(names)
    e = a.envelopes.(names{j});
    fprintf('%s %.0e %-9s %7.3f %.3f %.3e %.3f\n', ...
        label, a.sigma, names{j}, e.algebraic_slope, ...
        e.R2_algebraic, e.exponential_slope, e.R2_exponential);
end
end

function printMeshTable(records)
fprintf([' N     h       width  points lowerPts slopePts L2res ' ...
    'tailL2 totalL2 rate m_rho PG\n']);
for j = 1:numel(records)
    r = records(j);
    fprintf(['%4d %.3e %.3e %7.2f %8.2f %8.3g %.3e %.3e ' ...
        '%.3e %5.2f %6.2f %.2e\n'], ...
        r.N, r.h, r.width, r.points, r.censored_points_lower_bound, ...
        r.slope_points, r.resolved_L2, r.tail_L2, r.total_L2, ...
        r.state_rate, r.m_rho, r.final_pg);
end
end

function empirical = inferSigmaMinimum(fit, Nvalues, L, targetPoints)
empirical.N = Nvalues(:);
empirical.target_points = targetPoints;
empirical.sigma_min = NaN(size(empirical.N));
empirical.alpha = fit.alpha;
empirical.C = fit.C;
if isfinite(fit.alpha) && fit.alpha > 0 && isfinite(fit.C) && fit.C > 0
    h = 2 * L ./ empirical.N;
    empirical.sigma_min = (targetPoints * h / fit.C) .^ (1 / fit.alpha);
end
end

function printSigmaMinimum(empirical)
fprintf('\n------------------------------------------------------------\n');
fprintf(['Empirical minimum sigma for >=%g steepness-width points ' ...
    '(not a theorem)\n'], empirical.target_points);
fprintf('fit: width_slope = %.3e sigma^(%.6f)\n', ...
    empirical.C, empirical.alpha);
fprintf(' N       sigma_min\n');
for j = 1:numel(empirical.N)
    fprintf('%4d   %.3e\n', empirical.N(j), empirical.sigma_min(j));
end
end

function diagnosis = classifyLayerDiagnosis(fixed, mesh, scaling)
secondaryAlpha = scaling.secondaryH.alpha;
primaryCensored = all(cellfun(@(a) ~a.layer.primary.valid, fixed.analysisH));
ratesP2 = [mesh{1}.state_rate];
ratesP3 = [mesh{2}.state_rate];
ratesP4 = [mesh{3}.state_rate];
if any(ratesP2 >= 4) && ~any(ratesP3 >= 4) && ~any(ratesP4 >= 4)
    diagnosis.case = 'CASE L5';
    diagnosis.message = ['sigma=1e-6 enters a rapid-decay regime as the measured ' ...
        'steepness scale gains grid points, while sigma=1e-9/1e-12 remain in the ' ...
        'algebraic pre-asymptotic regime. Sigma is too small relative to ' ...
        'the tested spatial resolution. The prescribed q=0.1--0.9 width ' ...
        'is additionally box-censored and must not be used as a finite width.'];
    return;
end
if primaryCensored && isfinite(secondaryAlpha) && abs(secondaryAlpha) < 0.1
    diagnosis.case = 'CASE L4';
    diagnosis.message = ['The prescribed q=0.1--0.9 layer is box-censored ' ...
        'and the measurable q=0.25--0.75 width does not shrink materially ' ...
        'with sigma. Reject that full-band shrinking-width mechanism; the ' ...
        'steep inner scale remains a separate resolution diagnostic.'];
    return;
end
diagnosis.case = 'CASE L3';
diagnosis.message = ['Smooth-periodic-V state errors remain algebraic ' ...
    'without a resolved-band acceleration; the unresolved-layer ' ...
    'hypothesis is insufficient by itself.'];
end

function makeFigures(fixed, mesh, folder, sigmaList)
f1 = figure('Visible', 'off', 'Name', 'Sigma layer profiles');
tiledlayout(1, 2);
nexttile; hold on;
for index = 1:4
    a = fixed.analysisH{index};
    center = a.layer.primary.x_mid_right;
    mask = fixed.H{index}.grid.x >= 0;
    plot(fixed.H{index}.grid.x(mask) - center, ...
        a.layer.q_sigma(mask), 'DisplayName', ...
        sprintf('sigma=%.0e', sigmaList(index)));
end
xlim([-0.5, 2]); ylim([0, 1.05]); grid on;
xlabel('x-x_{q=0.5}'); ylabel('q_sigma'); title('harmonic V');
legend('Location', 'best');
nexttile; hold on;
for index = 2:4
    a = fixed.analysisP{index};
    center = a.layer.primary.x_mid_right;
    mask = fixed.P{index}.grid.x >= 0;
    plot(fixed.P{index}.grid.x(mask) - center, ...
        a.layer.q_sigma(mask), 'DisplayName', ...
        sprintf('sigma=%.0e', sigmaList(index)));
end
xlim([-0.5, 2]); ylim([0, 1.05]); grid on;
xlabel('x-x_{q=0.5}'); ylabel('q_sigma'); title('periodic V');
legend('Location', 'best');
exportgraphics(f1, fullfile(folder, 'sigma_layer_profiles.png'), ...
    'Resolution', 180);
close(f1);

f2 = figure('Visible', 'off', 'Name', 'Collapsed steep layer');
tiledlayout(1, 2);
nexttile; hold on;
for index = 1:4
    a = fixed.analysisH{index};
    center = a.layer.primary.x_mid_right;
    mask = fixed.H{index}.grid.x >= 0;
    coordinate = (fixed.H{index}.grid.x(mask) - center) ...
        / a.layer.width_slope;
    plot(coordinate, a.layer.q_sigma(mask), ...
        'DisplayName', sprintf('sigma=%.0e', sigmaList(index)));
end
xlim([-10, 20]); ylim([0, 1.05]); grid on;
xlabel('(x-x_{q=0.5})/width_{slope}'); ylabel('q_sigma');
title('harmonic V (primary width censored)'); legend('Location', 'best');
nexttile; hold on;
for index = 2:4
    a = fixed.analysisP{index};
    center = a.layer.primary.x_mid_right;
    mask = fixed.P{index}.grid.x >= 0;
    coordinate = (fixed.P{index}.grid.x(mask) - center) ...
        / a.layer.width_slope;
    plot(coordinate, a.layer.q_sigma(mask), ...
        'DisplayName', sprintf('sigma=%.0e', sigmaList(index)));
end
xlim([-10, 20]); ylim([0, 1.05]); grid on;
xlabel('(x-x_{q=0.5})/width_{slope}'); ylabel('q_sigma');
title('periodic V (primary width censored)'); legend('Location', 'best');
exportgraphics(f2, fullfile(folder, 'sigma_layer_collapsed.png'), ...
    'Resolution', 180);
close(f2);

f3 = figure('Visible', 'off', 'Name', 'Periodic sigma mesh errors');
hold on;
for j = 1:numel(mesh)
    loglog([mesh{j}.N], [mesh{j}.total_L2], 'o-', ...
        'DisplayName', sprintf('sigma=%.0e', mesh{j}(1).sigma));
end
grid on; xlabel('N'); ylabel('total spectral L2 error');
legend('Location', 'best');
exportgraphics(f3, fullfile(folder, 'sigma_layer_mesh_errors.png'), ...
    'Resolution', 180);
close(f3);

f4 = figure('Visible', 'off', 'Name', 'Layer widths versus sigma');
loglog(sigmaList, cellfun(@(a) a.layer.width_slope, fixed.analysisH), ...
    'o-', sigmaList, cellfun(@(a) a.layer.secondary.width_mean, ...
    fixed.analysisH), 's-');
grid on; xlabel('sigma'); ylabel('width diagnostic');
legend('steepness width', 'q=0.25--0.75 width', 'Location', 'best');
exportgraphics(f4, fullfile(folder, 'sigma_layer_width_scaling.png'), ...
    'Resolution', 180);
close(f4);
end
