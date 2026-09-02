%RUN_DIAG_FIXED_SIGMA_PERIODIC_MESH Fixed-sigma periodic-potential study.
%
% This is a diagnostic experiment. It changes neither the production
% harmonic potential nor any production solver or mesh-refinement code.
clearvars; clc;

% ====================== fixed mathematical problem =====================
epsilon = 1e-2;
sigma = 2e-6;

s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho)=rho+epsilon';

p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
potential_label = ...
    'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma, fixed sigma=2e-6';

beta = 10;
delta = 10;
mass = 1;
L = 8;
V_periodic = @(x) (L ^ 2 / pi ^ 2) .* (1 - cos(pi * x / L));
periodic_potential_label = ...
    'V_per(x)=(L^2/pi^2)*(1-cos(pi*x/L))';

N_list = [32, 64, 128, 256, 512, 1024, 2048];
N_ref = 4096;
% =======================================================================

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
reference_pg_tol = 1e-10;
reference_fourier_tol = 1e-12;
reference_tail_mass_tol = 1e-5;
reference_tail_max_tol = 1e-5;
reference_edge_density_tol = 1e-5;
reference_edge_derivative_tol = 1e-8;
reference_layer_points_tol = 8;
endpoint_jump_tol = 1e-10;
rate_floor = 1e-13;
save_result = true;
save_figures = true;
reuse_checkpoint = true;

allN = [N_list, N_ref];
if any(mod(allN, 2) ~= 0) || any(diff(allN) <= 0) ...
        || any(mod(N_ref, N_list) ~= 0) ...
        || any(mod(N_list(2:end), N_list(1:end-1)) ~= 0)
    error(['Use strictly increasing, nested, even Fourier grids with ' ...
        'N_ref divisible by every coarse N.']);
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
outputFolder = fullfile(root, 'results', ...
    'fixed_sigma_periodic_mesh');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
checkpointFile = fullfile(outputFolder, ...
    'fixed_sigma_2em6_periodic_checkpoint.mat');
outputFile = fullfile(outputFolder, ...
    'fixed_sigma_2em6_periodic_Nref4096.mat');

solver = makeSolver(main_pg_switch_tol, final_pg_tol);
problemTemplate.epsilon = epsilon;
problemTemplate.sigma = sigma;
problemTemplate.beta = beta;
problemTemplate.delta = delta;
problemTemplate.mass = mass;
problemTemplate.L = L;
problemTemplate.s_epsilon = s_epsilon;
problemTemplate.ds_epsilon = ds_epsilon;
problemTemplate.d2s_epsilon = d2s_epsilon;
problemTemplate.fisher_label = fisher_label;
problemTemplate.p_sigma = p_sigma;
problemTemplate.dp_sigma = dp_sigma;
problemTemplate.d2p_sigma = d2p_sigma;
problemTemplate.potential_label = potential_label;
problemTemplate.V_periodic = V_periodic;
problemTemplate.periodic_potential_label = periodic_potential_label;

solutions = cell(numel(allN), 1);
if reuse_checkpoint && isfile(checkpointFile)
    cached = load(checkpointFile, 'checkpoint');
    if isCompatibleCheckpoint(cached.checkpoint, allN, epsilon, sigma, L)
        solutions = cached.checkpoint.solutions;
        fprintf('Reused compatible fixed-sigma periodic checkpoint.\n');
    else
        warning('Ignoring incompatible fixed-sigma periodic checkpoint.');
    end
end

fprintf('\nFixed sigma = %.3e, smooth periodic potential\n', sigma);
fprintf('epsilon %.3e, beta %g, delta %g, L %g, Nref %d\n', ...
    epsilon, beta, delta, L, N_ref);

for index = 1:numel(allN)
    N = allN(index);
    if ~isempty(solutions{index})
        fprintf('N=%d reused (PG %.3e).\n', N, ...
            solutions{index}.diagnostics.final_pg_residual);
        continue;
    end
    grid = makeGrid(L, N);
    problem = makeProblem(grid, problemTemplate);
    if index == 1
        initial = loadCertifiedFixedSigmaWarmStart( ...
            root, grid, epsilon, sigma, mass);
        if isempty(initial)
            initial = mass * ones(grid.N, 1) / grid.domain_length;
            fprintf(['  no certified fixed-sigma warm state found; ' ...
                'using a strictly positive constant state.\n']);
        end
    else
        previous = solutions{index - 1};
        initial = src.discretization.ps.Prolong(previous.rho, N);
        initial = src.constraints.ProjectPositiveConservative( ...
            initial, mass, grid.h, solver.projection_tol);
    end
    fprintf('Solving N=%d ...\n', N);
    timer = tic;
    result = src.SolveGroundState1D(problem, initial, solver);
    result.grid = grid;
    result.solve_wall_time = toc(timer);
    solutions{index} = result;
    fprintf(['  done in %.2f s: E %.15e, PG %.3e, KKT %.3e, ' ...
        'main/polish %d/%d\n'], result.solve_wall_time, ...
        result.target_energy, result.diagnostics.final_pg_residual, ...
        result.diagnostics.final_kkt_residual, ...
        result.diagnostics.main_iterations, ...
        result.diagnostics.polish_iterations);

    checkpoint.N_values = allN;
    checkpoint.epsilon = epsilon;
    checkpoint.sigma = sigma;
    checkpoint.L = L;
    checkpoint.solutions = solutions;
    save(checkpointFile, 'checkpoint', '-v7.3');
end

analysis = cell(numel(allN), 1);
for index = 1:numel(allN)
    analysis{index} = analyzeState(solutions{index});
end
reference = solutions{end};
referenceAnalysis = analysis{end};

optimizationPass = reference.diagnostics.final_pg_residual ...
    <= reference_pg_tol && reference.diagnostics.exact_zero_count == 0;
fourierPass = referenceAnalysis.fourier.tail_ratio_quarter ...
    <= reference_fourier_tol;
boxPass = referenceAnalysis.far_field.tail_mass ...
    <= reference_tail_mass_tol ...
    && referenceAnalysis.far_field.tail_max <= reference_tail_max_tol ...
    && referenceAnalysis.far_field.edge_density ...
    <= reference_edge_density_tol ...
    && referenceAnalysis.far_field.edge_abs_drho ...
    <= reference_edge_derivative_tol;
layerPass = referenceAnalysis.layer.points_per_slope_width ...
    >= reference_layer_points_tol;
endpointPass = abs(referenceAnalysis.endpoint.GV_value_jump) ...
    <= endpoint_jump_tol ...
    && abs(referenceAnalysis.endpoint.GV_derivative_jump) ...
    <= endpoint_jump_tol;
referenceFullyCertified = optimizationPass && fourierPass ...
    && boxPass && layerPass && endpointPass;

fprintf('\nReference diagnostics (Nref=%d)\n', N_ref);
fprintf('  final PG           : %.3e\n', ...
    reference.diagnostics.final_pg_residual);
fprintf('  KKT residual       : %.3e\n', ...
    reference.diagnostics.final_kkt_residual);
fprintf('  Fourier quarter    : %.3e\n', ...
    referenceAnalysis.fourier.tail_ratio_quarter);
fprintf('  tail mass/max      : %.3e / %.3e\n', ...
    referenceAnalysis.far_field.tail_mass, ...
    referenceAnalysis.far_field.tail_max);
fprintf('  edge rho/|Drho|    : %.3e / %.3e\n', ...
    referenceAnalysis.far_field.edge_density, ...
    referenceAnalysis.far_field.edge_abs_drho);
fprintf('  exact zero count   : %d\n', ...
    reference.diagnostics.exact_zero_count);
fprintf('  steep-width points : %.3f\n', ...
    referenceAnalysis.layer.points_per_slope_width);
fprintf('  optimization reference : %s\n', passFail(optimizationPass));
fprintf('  Fourier reference      : %s\n', passFail(fourierPass));
fprintf('  box reference          : %s\n', passFail(boxPass));
fprintf('  layer resolution       : %s\n', passFail(layerPass));
fprintf('  periodic GV endpoint   : %s\n', passFail(endpointPass));
if ~referenceFullyCertified
    warning('N_ref=4096 is not fully certified.');
end

records = repmat(struct(), numel(N_list), 1);
for index = 1:numel(N_list)
    comparison = src.diagnostics.SpectralStateComparison( ...
        solutions{index}.rho, solutions{index}.grid, ...
        reference.rho, reference.grid);
    records(index).N = N_list(index);
    records(index).h = solutions{index}.grid.h;
    records(index).steepness_width = ...
        analysis{index}.layer.width_slope;
    records(index).points_per_steepness_width = ...
        analysis{index}.layer.points_per_slope_width;
    records(index).max_abs_dqdx = analysis{index}.layer.max_abs_dqdx;
    records(index).full_crossing_valid = ...
        analysis{index}.layer.primary.valid;
    records(index).full_crossing_width = ...
        analysis{index}.layer.primary.width_mean;
    records(index).resolved_L2 = comparison.resolved_L2_error;
    records(index).resolved_Linf = comparison.resolved_Linf_error;
    records(index).tail_L2 = comparison.reference_tail_L2;
    records(index).total_L2 = comparison.total_spectral_L2_error;
    records(index).target_energy_error = abs( ...
        solutions{index}.target_energy - reference.target_energy);
    records(index).baseline_energy_difference = ...
        solutions{index}.baseline_energy - reference.baseline_energy;
    records(index).final_pg = ...
        solutions{index}.diagnostics.final_pg_residual;
    records(index).final_kkt = ...
        solutions{index}.diagnostics.final_kkt_residual;
    records(index).fft_tail = ...
        analysis{index}.fourier.tail_ratio_quarter;
    records(index).envelope = analysis{index}.envelope;
    records(index).envelope_regime = analysis{index}.envelope_regime;
    records(index).endpoint = analysis{index}.endpoint;
end

ratesResolved = rateSeries([records.resolved_L2], rate_floor);
ratesTail = rateSeries([records.tail_L2], rate_floor);
ratesTotal = rateSeries([records.total_L2], rate_floor);
ratesEnergy = rateSeries([records.target_energy_error], rate_floor);
for index = 1:numel(records)
    records(index).rate_resolved = ratesResolved(index);
    records(index).rate_tail = ratesTail(index);
    records(index).rate_total = ratesTotal(index);
    records(index).rate_energy = ratesEnergy(index);
end

printLayerTable(allN, solutions, analysis);
printErrorTables(records, rate_floor);
printCorrelationTable(records, rate_floor);
printEnvelopeTable(allN, analysis);
printEnvelopeBands(allN, analysis);
printEndpointTable(allN, analysis);

rapidAcceleration = detectRapidAcceleration(records, analysis);
envelopeTransition = detectEnvelopeTransition(analysis);
diagnosis = classifyDiagnosis(rapidAcceleration, envelopeTransition, ...
    referenceFullyCertified, layerPass, endpointPass);
fprintf('\n============================================================\n');
fprintf('DIAGNOSIS %s\n%s\n', diagnosis.case, diagnosis.message);
fprintf('rapid acceleration : %s\n', passFail(rapidAcceleration));
fprintf('envelope transition: %s\n', passFail(envelopeTransition));
fprintf('============================================================\n');

archive.description = ['Diagnostic fixed-sigma smooth-periodic-potential ' ...
    'mesh refinement; not a production harmonic-potential result.'];
archive.epsilon = epsilon;
archive.sigma = sigma;
archive.beta = beta;
archive.delta = delta;
archive.mass = mass;
archive.L = L;
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
archive.periodic_potential.V = V_periodic;
archive.periodic_potential.label = periodic_potential_label;
archive.N_list = N_list;
archive.N_ref = N_ref;
archive.grid = cellfun(@(state) state.grid, solutions, ...
    'UniformOutput', false);
archive.rho = cellfun(@(state) state.rho, solutions, ...
    'UniformOutput', false);
archive.target_energy = cellfun(@(state) state.target_energy, solutions);
archive.baseline_energy = cellfun(@(state) state.baseline_energy, solutions);
archive.records = records;
archive.analysis = analysis;
archive.solver_diagnostics = cellfun(@(state) state.diagnostics, ...
    solutions, 'UniformOutput', false);
archive.solver = solver;
archive.reference.optimization_pass = optimizationPass;
archive.reference.fourier_pass = fourierPass;
archive.reference.box_pass = boxPass;
archive.reference.layer_pass = layerPass;
archive.reference.endpoint_pass = endpointPass;
archive.reference.fully_certified = referenceFullyCertified;
archive.rapid_acceleration = rapidAcceleration;
archive.envelope_transition = envelopeTransition;
archive.diagnosis = diagnosis;
if save_result
    save(outputFile, 'archive', '-v7.3');
end
if save_figures
    makeFigures(records, allN, analysis, outputFolder, rate_floor);
end

function solver = makeSolver(mainTolerance, finalTolerance)
config = experiments.DefaultConfig();
solver = config.solver;
solver.name = 'fista_cd';
solver.splitting = 'potential_prox';
solver.projection_name = 'semismooth';
solver.projection_tol = 1e-14;
solver.pg_tol = mainTolerance;
solver.final_pg_tol = finalTolerance;
solver.certification_tol = finalTolerance;
solver.residual_check_interval = 10;
solver.potential_prox.mass_tol = 1e-14;
solver.potential_prox.inner_tol = 1e-14;
% Early coarse-grid candidates can create a very wide safeguarded node
% bracket. This changes only the iteration safety cap, not the prox
% equation or either tolerance.
solver.potential_prox.inner_max_iter = 100;
solver.max_iter = 200000;
solver.switch.enabled = true;
solver.switch.energy_window = 50;
solver.switch.energy_tol = 1e-12;
solver.switch.consecutive_windows = 2;
solver.switch.min_iter = 200;
solver.switch.pg_entry_tol = mainTolerance;
solver.switch.max_main_iter = 20000;
solver.polish_mode = 'if_needed';
solver.polish.pg_tol = finalTolerance;
solver.polish.max_iter = 20;
solver.polish.linear_solver = 'interior_pcg_schur';
solver.polish.preconditioner = 'fd_variable';
solver.polish.allow_pdas_fallback = true;
solver.display = false;
solver.display_every = 200;
end

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end

function problem = makeProblem(grid, template)
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = template.V_periodic(grid.x);
problem.beta = template.beta;
problem.delta = template.delta;
problem.mass = template.mass;
problem.fisher_regularization.epsilon = template.epsilon;
problem.fisher_regularization.s_epsilon = template.s_epsilon;
problem.fisher_regularization.ds_epsilon = template.ds_epsilon;
problem.fisher_regularization.d2s_epsilon = template.d2s_epsilon;
problem.fisher_regularization.label = template.fisher_label;
problem.potential_regularization.name = 'inline_fixed_sigma_periodic_diag';
problem.potential_regularization.sigma = template.sigma;
problem.potential_regularization.p_sigma = template.p_sigma;
problem.potential_regularization.dp_sigma = template.dp_sigma;
problem.potential_regularization.d2p_sigma = template.d2p_sigma;
problem.potential_regularization.label = template.potential_label;
problem.potential_regularization.prox_type = 'generic_convex';
problem.potential_regularization.convexity_tol = 1e-14;
problem.potential_label = template.periodic_potential_label;
problem.entropy.enabled = false;
problem.entropy.eta = 0;
end

function initial = loadCertifiedFixedSigmaWarmStart( ...
    root, grid, epsilon, sigma, mass)
initial = [];
filename = fullfile(root, 'results', 'potential_mesh', ...
    'potential_mesh_inline_p_sigma_sigma_2em06_Nref8192_eps_1em03.mat');
if ~isfile(filename) || epsilon ~= 1e-3 || sigma ~= 2e-6
    return;
end
data = load(filename, 'rho', 'grid');
Nvalues = cellfun(@(item) item.N, data.grid);
index = find(Nvalues == grid.N, 1);
if isempty(index)
    return;
end
candidate = data.rho{index}(:);
if any(~isfinite(candidate)) || any(candidate <= 0)
    return;
end
initial = candidate * (mass / src.constraints.Mass(candidate, grid.h));
fprintf(['  using certified same-(epsilon,sigma,N) positive state as ' ...
    'the first periodic-objective warm start.\n']);
end

function compatible = isCompatibleCheckpoint(checkpoint, Nvalues, ...
    epsilon, sigma, L)
required = {'N_values', 'epsilon', 'sigma', 'L', 'solutions'};
compatible = isstruct(checkpoint) && all(isfield(checkpoint, required)) ...
    && isequal(checkpoint.N_values, Nvalues) ...
    && checkpoint.epsilon == epsilon && checkpoint.sigma == sigma ...
    && checkpoint.L == L && numel(checkpoint.solutions) == numel(Nvalues);
end

function diagnostic = analyzeState(result)
diagnostic.layer = src.diagnostics.TransitionLayerDiagnostics( ...
    result.rho, result.grid, result.problem.plan, ...
    result.problem.potential_regularization);
diagnostic.fourier = src.diagnostics.FourierTailDiagnostics(result.rho);
diagnostic.far_field = src.diagnostics.FarFieldDiagnostics( ...
    result.rho, result.grid, result.problem.plan);
fitOptions.bands = [8, 16; 16, 32; 32, 64; 64, 128; ...
    128, 256; 256, 512; 512, 1024; 1024, 2048];
fitOptions.fit_range = [8, min(2048, result.grid.N / 2 - 1)];
fitOptions.relative_floor = 100 * eps;
diagnostic.envelope = src.diagnostics.FourierEnvelopeFit( ...
    result.rho, fitOptions);
diagnostic.envelope_regime = envelopeRegime(diagnostic.envelope);
L = result.grid.L;
endpointData.V = [2 * L ^ 2 / pi ^ 2, 2 * L ^ 2 / pi ^ 2];
endpointData.dV = [0, 0];
diagnostic.endpoint = src.diagnostics.ActualEndpointCompatibility( ...
    result.rho, result.problem, endpointData);
end

function label = envelopeRegime(envelope)
validBands = envelope.band_max(isfinite(envelope.band_max));
if envelope.fit_count < 4 || (numel(validBands) >= 2 ...
        && all(validBands(max(1, end - 1):end) ...
        <= 1e-13 * max(1, max(envelope.amplitudes))))
    label = 'roundoff/floor';
elseif envelope.R2_exponential >= envelope.R2_algebraic + 0.03
    label = 'exponential-looking';
elseif envelope.R2_algebraic >= envelope.R2_exponential + 0.03
    label = 'algebraic-looking';
else
    label = 'mixed/pre-asymptotic';
end
end

function rates = rateSeries(errors, floorValue)
errors = errors(:).';
rates = NaN(size(errors));
for index = 2:numel(errors)
    if isfinite(errors(index - 1)) && isfinite(errors(index)) ...
            && errors(index - 1) > floorValue && errors(index) > floorValue ...
            && errors(index) > 0
        rates(index) = log(errors(index - 1) / errors(index)) / log(2);
    end
end
end

function printLayerTable(Nvalues, solutions, analysis)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Table A: layer resolution and solver certification\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' N       h          steep_width   pts_layer  max|dqdx|  PG_final  FFT_tail  label\n');
for index = 1:numel(Nvalues)
    layer = analysis{index}.layer;
    fprintf('%4d  %.3e  %.3e  %9.3f  %.3e  %.2e  %.3e  %s\n', ...
        Nvalues(index), solutions{index}.grid.h, layer.width_slope, ...
        layer.points_per_slope_width, layer.max_abs_dqdx, ...
        solutions{index}.diagnostics.final_pg_residual, ...
        analysis{index}.fourier.tail_ratio_quarter, ...
        layer.slope_resolution_label);
end
end

function printErrorTables(records, floorValue)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Table B: spectral state and target-energy errors\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' N      L2_res       rate       L2_tail      rate       L2_total     rate       dE           rate\n');
for index = 1:numel(records)
    r = records(index);
    fprintf('%4d  %.3e  %8s  %.3e  %8s  %.3e  %8s  %.3e  %8s\n', ...
        r.N, r.resolved_L2, rateLabel(r.rate_resolved, ...
        r.resolved_L2, floorValue), r.tail_L2, ...
        rateLabel(r.rate_tail, r.tail_L2, floorValue), ...
        r.total_L2, rateLabel(r.rate_total, r.total_L2, floorValue), ...
        r.target_energy_error, rateLabel(r.rate_energy, ...
        r.target_energy_error, floorValue));
end
end

function printCorrelationTable(records, floorValue)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Layer-resolution / state-convergence correlation\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' N     points/layer    L2_total       state_rate\n');
for index = 1:numel(records)
    r = records(index);
    fprintf('%4d    %10.3f    %.3e    %s\n', r.N, ...
        r.points_per_steepness_width, r.total_L2, ...
        rateLabel(r.rate_total, r.total_L2, floorValue));
end
end

function printEnvelopeTable(Nvalues, analysis)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Fourier envelope fits\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' N      m_alg   R2_alg    c_exp      R2_exp   fit_count  regime\n');
for index = 1:numel(Nvalues)
    envelope = analysis{index}.envelope;
    fprintf('%4d  %8.3f  %.4f  %.3e  %.4f  %9d  %s\n', ...
        Nvalues(index), envelope.algebraic_slope, ...
        envelope.R2_algebraic, envelope.exponential_slope, ...
        envelope.R2_exponential, envelope.fit_count, ...
        analysis{index}.envelope_regime);
end
end

function printEnvelopeBands(Nvalues, analysis)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Fourier coefficient band envelopes\n');
fprintf('----------------------------------------------------------------------\n');
for index = 1:numel(Nvalues)
    envelope = analysis{index}.envelope;
    fprintf('N=%d\n', Nvalues(index));
    for band = 1:size(envelope.bands, 1)
        if envelope.band_sample_count(band) > 0 ...
                && isfinite(envelope.band_max(band))
            fprintf('  |k|=%4d:%-4d  max %.3e  rms %.3e\n', ...
                envelope.bands(band, 1), envelope.bands(band, 2), ...
                envelope.band_max(band), envelope.band_rms(band));
        end
    end
end
end

function printEndpointTable(Nvalues, analysis)
fprintf('\n----------------------------------------------------------------------\n');
fprintf('Smooth-periodic Euler-potential endpoint compatibility\n');
fprintf('----------------------------------------------------------------------\n');
fprintf(' N      GV_value_jump     GV_derivative_jump\n');
for index = 1:numel(Nvalues)
    endpoint = analysis{index}.endpoint;
    fprintf('%4d   %+.3e          %+.3e\n', Nvalues(index), ...
        endpoint.GV_value_jump, endpoint.GV_derivative_jump);
end
end

function label = rateLabel(rate, errorValue, floorValue)
if errorValue <= floorValue
    label = 'floor';
elseif isfinite(rate)
    label = sprintf('%.2f', rate);
else
    label = '--';
end
end

function accelerated = detectRapidAcceleration(records, analysis)
rates = [records.rate_total];
finiteRates = rates(isfinite(rates));
points = cellfun(@(item) item.layer.points_per_slope_width, analysis);
if numel(finiteRates) < 3
    accelerated = false;
    return;
end
rateGain = max(finiteRates) - median(finiteRates(1:min(2, end)));
growthCount = nnz(diff(finiteRates) > 0.5);
crossesLayerRange = min(points) < 4 && max(points) >= 8;
accelerated = max(finiteRates) >= 4 && rateGain >= 2 ...
    && growthCount >= 2 && crossesLayerRange;
end

function transitioned = detectEnvelopeTransition(analysis)
labels = cellfun(@(item) item.envelope_regime, analysis, ...
    'UniformOutput', false);
earlyAlgebraic = any(strcmp(labels(1:min(4, end)), ...
    'algebraic-looking')) || any(strcmp(labels(1:min(4, end)), ...
    'mixed/pre-asymptotic'));
lateRapid = any(strcmp(labels(max(1, end - 2):end), ...
    'exponential-looking')) || any(strcmp(labels(max(1, end - 2):end), ...
    'roundoff/floor'));
transitioned = earlyAlgebraic && lateRapid;
end

function diagnosis = classifyDiagnosis(acceleration, envelopeTransition, ...
    referenceCertified, layerPass, endpointPass)
if ~endpointPass
    diagnosis.case = 'CASE S4';
    diagnosis.message = ['Periodic GV endpoint compatibility failed; ' ...
        'stop spatial-convergence interpretation.'];
elseif acceleration && envelopeTransition && referenceCertified
    diagnosis.case = 'CASE S1';
    diagnosis.message = ['Once the fixed-sigma steep layer is resolved, ' ...
        'state errors accelerate and the Fourier envelope reaches an ' ...
        'exponential/floor regime.'];
elseif acceleration
    diagnosis.case = 'CASE S3';
    diagnosis.message = ['Evidence supports spectral acceleration, but ' ...
        'N_ref=4096 is not fully certified for quantitative confirmation.'];
elseif layerPass
    diagnosis.case = 'CASE S2';
    diagnosis.message = ['The reference layer has at least eight points, ' ...
        'but state convergence does not show sustained acceleration; the ' ...
        'unresolved sigma layer is not a sufficient explanation.'];
else
    diagnosis.case = 'CASE S3';
    diagnosis.message = ['N_ref=4096 does not resolve the fixed-sigma ' ...
        'layer sufficiently, so the proposed acceleration cannot yet be ' ...
        'confirmed or rejected.'];
end
end

function makeFigures(records, Nvalues, analysis, folder, floorValue)
Ncoarse = [records.N];
f1 = figure('Visible', 'off', 'Name', 'Fixed sigma state errors');
loglog(Ncoarse, max([records.resolved_L2], floorValue), 'o-', ...
    Ncoarse, max([records.tail_L2], floorValue), 's-', ...
    Ncoarse, max([records.total_L2], floorValue), 'd-');
grid on; xlabel('N'); ylabel('L2 error');
legend('resolved L2', 'reference tail L2', 'total spectral L2', ...
    'Location', 'best');
exportgraphics(f1, fullfile(folder, '01_state_errors.png'), ...
    'Resolution', 180);
close(f1);

f2 = figure('Visible', 'off', 'Name', 'Fixed sigma energy error');
loglog(Ncoarse, max([records.target_energy_error], floorValue), 'o-');
grid on; xlabel('N'); ylabel('|E_N-E_{ref}|');
exportgraphics(f2, fullfile(folder, '02_target_energy_error.png'), ...
    'Resolution', 180);
close(f2);

f3 = figure('Visible', 'off', 'Name', 'Fixed sigma layer resolution');
points = cellfun(@(item) item.layer.points_per_slope_width, analysis);
semilogx(Nvalues, points, 'o-'); hold on;
yline(4, '--', '4 points'); yline(8, '--', '8 points');
grid on; xlabel('N'); ylabel('points per steepness width');
exportgraphics(f3, fullfile(folder, '03_layer_resolution.png'), ...
    'Resolution', 180);
close(f3);

representative = [128, 512, 1024, 2048, 4096];
f4 = figure('Visible', 'off', 'Name', 'Fixed sigma Fourier envelopes');
hold on;
for value = representative
    index = find(Nvalues == value, 1);
    envelope = analysis{index}.envelope;
    [absoluteModes, order] = sort(abs(envelope.modes));
    semilogy(absoluteModes, max(envelope.amplitudes(order), realmin), ...
        '.', 'DisplayName', sprintf('N=%d', value));
end
grid on; xlabel('|k|'); ylabel('|rho hat_k|');
legend('Location', 'best');
exportgraphics(f4, fullfile(folder, '04_fourier_envelopes.png'), ...
    'Resolution', 180);
close(f4);
end

function label = passFail(value)
if value
    label = 'PASS';
else
    label = 'FAIL';
end
end
