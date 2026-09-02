%RUN_MESH_REFINEMENT_SPECTRAL_HARMONIC Formal Fourier mesh experiment.
%
% The physical trap remains exactly harmonic on |x|<=R0. Only the
% negligible-density boundary strip is continued to a high constant using
% a flat C-infinity step, making the periodic extension smooth. This run
% does not alter the production harmonic solver or any other experiment.
clearvars; clc;

% ====================== fixed mathematical problem =====================
epsilon = 1e-2;
sigma = 2e-4;
beta = 10;
delta = 10;
mass = 1;
L = 8;
R0 = 0.75 * L;
R1 = 0.90 * L;

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

N_list = [32, 64, 128, 256, 512];
N_ref = 1024;
% =======================================================================

main_pg_switch_tol = 1e-5;
final_pg_tol = 1e-12;
rate_floor = 1e-13;
transition_density_warning_tol = 1e-6;
save_result = true;
save_figures = true;

allN = [N_list, N_ref];
if any(mod(allN, 2) ~= 0) || any(diff(allN) <= 0) ...
        || any(mod(N_ref, N_list) ~= 0) ...
        || any(mod(N_list(2:end), N_list(1:end-1)) ~= 0)
    error('Use strictly nested even Fourier grids.');
end
if ~(0 < R0 && R0 < R1 && R1 < L)
    error('Require 0<R0<R1<L.');
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
outputFolder = fullfile(root, 'results', ...
    'mesh_refinement_spectral_harmonic');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
outputFile = fullfile(outputFolder, ...
    ['mesh_refinement_harmonic_periodicized_' ...
    'eps1em2_sigma2em4_Nref1024.mat']);

solver = makeSolver(main_pg_switch_tol, final_pg_tol);
template.epsilon = epsilon;
template.sigma = sigma;
template.beta = beta;
template.delta = delta;
template.mass = mass;
template.L = L;
template.R0 = R0;
template.R1 = R1;
template.s_epsilon = s_epsilon;
template.ds_epsilon = ds_epsilon;
template.d2s_epsilon = d2s_epsilon;
template.fisher_label = fisher_label;
template.p_sigma = p_sigma;
template.dp_sigma = dp_sigma;
template.d2p_sigma = d2p_sigma;
template.potential_label = potential_label;

solutions = cell(numel(allN), 1);
warmStart = [];
for index = 1:numel(allN)
    N = allN(index);
    grid = makeGrid(L, N);
    [V, potentialData] = periodicizedHarmonic(grid.x, L, R0, R1);
    validatePotential(V, potentialData, grid, R0);
    problem = makeProblem(grid, V, template, ...
        'harmonic trap with C-infinity boundary periodicization');
    if isempty(warmStart)
        warmStart = mass * ones(N, 1) / grid.domain_length;
    else
        warmStart = src.discretization.ps.Prolong(warmStart, N);
        warmStart = src.constraints.ProjectPositiveConservative( ...
            warmStart, mass, grid.h, solver.projection_tol);
    end
    fprintf('Solving periodicized harmonic N=%d ...\n', N);
    timer = tic;
    result = src.SolveGroundState1D(problem, warmStart, solver);
    result.grid = grid;
    result.solve_wall_time = toc(timer);
    result.fourier = src.diagnostics.FourierTailDiagnostics(result.rho);
    result.far_field = src.diagnostics.FarFieldDiagnostics( ...
        result.rho, grid, result.problem.plan);
    result.layer = src.diagnostics.TransitionLayerDiagnostics( ...
        result.rho, grid, result.problem.plan, ...
        result.problem.potential_regularization);
    solutions{index} = result;
    warmStart = result.rho;
    fprintf('  %.2f s, E %.15e, PG %.3e, KKT %.3e\n', ...
        result.solve_wall_time, result.target_energy, ...
        result.diagnostics.final_pg_residual, ...
        result.diagnostics.final_kkt_residual);
end

reference = solutions{end};
referenceGrid = reference.grid;
[Vreference, potentialValidation] = periodicizedHarmonic( ...
    referenceGrid.x, L, R0, R1);
Hreference = 0.5 * referenceGrid.x .^ 2;
dVreference = src.discretization.ps.FirstDerivative( ...
    Vreference, reference.problem.plan);
periodicValueJump = abs(Vreference(end) - Vreference(1));
periodicDerivativeJump = abs(dVreference(end) - dVreference(1));
maxInteriorDifference = max(abs( ...
    Vreference(abs(referenceGrid.x) <= R0) ...
    - Hreference(abs(referenceGrid.x) <= R0)));
maxPotentialModification = max(abs(Vreference - Hreference));
minPotential = min(Vreference);
if maxInteriorDifference > eps(max(1, max(abs(Hreference))))
    error('Periodicized potential changed the harmonic interior.');
end
if minPotential < 0
    error('Periodicized harmonic potential must remain nonnegative.');
end

transitionMask = abs(referenceGrid.x) >= R0;
rhoTransitionMax = max(reference.rho(transitionMask));
transitionMass = referenceGrid.h * sum(reference.rho(transitionMask));
deltaV = Vreference - Hreference;
potentialModificationEnergy = referenceGrid.h * sum( ...
    deltaV .* p_sigma(reference.rho));
absPotentialModificationEnergy = abs(potentialModificationEnergy);
if rhoTransitionMax > transition_density_warning_tol
    warning(['Harmonic potential is being modified in a ' ...
        'non-negligible-density region; increase L or R0.']);
end

records = repmat(struct(), numel(N_list), 1);
for index = 1:numel(N_list)
    comparison = src.diagnostics.SpectralStateComparison( ...
        solutions{index}.rho, solutions{index}.grid, ...
        reference.rho, reference.grid);
    records(index).N = N_list(index);
    records(index).resolved_L2 = comparison.resolved_L2_error;
    records(index).tail_L2 = comparison.reference_tail_L2;
    records(index).total_L2 = comparison.total_spectral_L2_error;
    records(index).resolved_Linf = comparison.resolved_Linf_error;
    records(index).target_energy_error = abs( ...
        solutions{index}.target_energy - reference.target_energy);
    records(index).final_pg = ...
        solutions{index}.diagnostics.final_pg_residual;
end
ratesL2 = rateSeries([records.total_L2], rate_floor);
ratesLinf = rateSeries([records.resolved_Linf], rate_floor);
ratesEnergy = rateSeries([records.target_energy_error], rate_floor);
for index = 1:numel(records)
    records(index).rate_L2 = ratesL2(index);
    records(index).rate_Linf = ratesLinf(index);
    records(index).rate_energy = ratesEnergy(index);
end

% One reference-grid sanity solve under the direct harmonic trap.
harmonicProblem = makeProblem(referenceGrid, Hreference, template, ...
    'direct harmonic diagnostic');
fprintf('Solving direct harmonic sanity state at Nref=%d ...\n', N_ref);
harmonicTimer = tic;
harmonicResult = src.SolveGroundState1D( ...
    harmonicProblem, reference.rho, solver);
harmonicResult.grid = referenceGrid;
harmonicResult.solve_wall_time = toc(harmonicTimer);
stateDifference = reference.rho - harmonicResult.rho;
globalHarmonicDifference = sqrt(referenceGrid.h ...
    * sum(stateDifference .^ 2));
centralMask = abs(referenceGrid.x) <= R0;
centralHarmonicDifference = sqrt(referenceGrid.h ...
    * sum(stateDifference(centralMask) .^ 2));
harmonicEnergyPeriodicizedState = src.discretization.ps.Energy( ...
    reference.rho, harmonicProblem);
harmonicEnergyDifference = abs( ...
    harmonicEnergyPeriodicizedState - harmonicResult.target_energy);

lastStateError = records(end).total_L2;
lastEnergyError = records(end).target_energy_error;
transitionMassRatio = transitionMass / max(lastStateError, realmin);
modificationEnergyRatio = absPotentialModificationEnergy ...
    / max(lastEnergyError, realmin);
boundaryModificationWarning = centralHarmonicDifference >= lastStateError;
if boundaryModificationWarning
    warning(['Boundary periodicization is not negligible: its central ' ...
        'state difference is at least the N=512 spatial error.']);
end

clearAcceleration = detectAcceleration(records);
printMainTable(records, epsilon, sigma, L, R0, R1, N_ref, rate_floor);
fprintf('\nReference summary\n');
fprintf('  PG final                       : %.3e\n', ...
    reference.diagnostics.final_pg_residual);
fprintf('  KKT final                      : %.3e\n', ...
    reference.diagnostics.final_kkt_residual);
fprintf('  Fourier tail                   : %.3e\n', ...
    reference.fourier.tail_ratio_quarter);
fprintf('  exact zeros                    : %d\n', ...
    reference.diagnostics.exact_zero_count);
fprintf('  points per sigma layer         : %.3f\n', ...
    reference.layer.points_per_slope_width);
fprintf('  periodic V value jump          : %.3e\n', periodicValueJump);
fprintf('  periodic V derivative jump     : %.3e\n', periodicDerivativeJump);
fprintf('  rho transition max             : %.3e\n', rhoTransitionMax);
fprintf('  transition mass                : %.3e\n', transitionMass);
fprintf('  potential modification energy  : %+.3e\n', ...
    potentialModificationEnergy);
fprintf('  max interior V difference      : %.3e\n', ...
    maxInteriorDifference);
fprintf('  min(V), max modification       : %.3e, %.3e\n', ...
    minPotential, maxPotentialModification);
fprintf('\nBoundary-modification diagnostics\n');
fprintf('  transition mass / N512 L2      : %.3e\n', ...
    transitionMassRatio);
fprintf('  |dE potential| / N512 dE       : %.3e\n', ...
    modificationEnergyRatio);
fprintf('  perH/direct-H global L2         : %.3e\n', ...
    globalHarmonicDifference);
fprintf('  perH/direct-H central L2        : %.3e\n', ...
    centralHarmonicDifference);
fprintf('  direct-H energy difference      : %.3e\n', ...
    harmonicEnergyDifference);
fprintf('  clear spectral-type acceleration: %s\n', ...
    passFail(clearAcceleration));

archive.description = ['Fourier mesh refinement for a harmonic trap ' ...
    'with C-infinity boundary-only periodic continuation.'];
archive.epsilon = epsilon;
archive.sigma = sigma;
archive.beta = beta;
archive.delta = delta;
archive.mass = mass;
archive.L = L;
archive.R0 = R0;
archive.R1 = R1;
archive.N_list = N_list;
archive.N_ref = N_ref;
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
archive.potential.V_periodicized_harmonic = Vreference;
archive.potential.V_harmonic = Hreference;
archive.potential.chi = potentialValidation.chi;
archive.potential.Vflat = potentialValidation.Vflat;
archive.potential.max_interior_difference = maxInteriorDifference;
archive.potential.max_modification = maxPotentialModification;
archive.potential.min_value = minPotential;
archive.grid = cellfun(@(state) state.grid, solutions, ...
    'UniformOutput', false);
archive.rho = cellfun(@(state) state.rho, solutions, ...
    'UniformOutput', false);
archive.target_energy = cellfun(@(state) state.target_energy, solutions);
archive.records = records;
archive.reference.diagnostics = reference.diagnostics;
archive.reference.fourier = reference.fourier;
archive.reference.layer = reference.layer;
archive.reference.rho_transition_max = rhoTransitionMax;
archive.reference.transition_mass = transitionMass;
archive.reference.potential_modification_energy = ...
    potentialModificationEnergy;
archive.reference.periodic_value_jump = periodicValueJump;
archive.reference.periodic_derivative_jump = periodicDerivativeJump;
archive.harmonic_sanity.rho = harmonicResult.rho;
archive.harmonic_sanity.diagnostics = harmonicResult.diagnostics;
archive.harmonic_sanity.global_L2_difference = globalHarmonicDifference;
archive.harmonic_sanity.central_L2_difference = centralHarmonicDifference;
archive.harmonic_sanity.harmonic_energy_difference = ...
    harmonicEnergyDifference;
archive.boundary_diagnostics.transition_mass_ratio = transitionMassRatio;
archive.boundary_diagnostics.modification_energy_ratio = ...
    modificationEnergyRatio;
archive.boundary_diagnostics.warning = boundaryModificationWarning;
archive.clear_spectral_acceleration = clearAcceleration;
archive.solver = solver;
if save_result
    save(outputFile, 'archive', '-v7.3');
end
writeLatexTable(records, fullfile(outputFolder, ...
    'mesh_refinement_table.tex'), rate_floor);
if save_figures
    makeFigures(solutions, records, Vreference, Hreference, ...
        referenceGrid, R0, R1, outputFolder, rate_floor);
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

function [V, diagnostic] = periodicizedHarmonic(x, L, R0, R1)
harmonic = 0.5 * x .^ 2;
t = (abs(x) - R0) / (R1 - R0);
chi = model.SmoothStepCInf(t);
Vflat = 0.5 * L ^ 2;
V = (1 - chi) .* harmonic + chi .* Vflat;
diagnostic.harmonic = harmonic;
diagnostic.chi = chi;
diagnostic.Vflat = Vflat;
end

function validatePotential(V, diagnostic, grid, R0)
interior = abs(grid.x) <= R0;
interiorDifference = max(abs(V(interior) ...
    - diagnostic.harmonic(interior)));
if interiorDifference > eps(max(1, max(abs(diagnostic.harmonic))))
    error('The periodicized potential changed V=x^2/2 inside |x|<=R0.');
end
if min(V) < 0 || any(~isfinite(V))
    error('The periodicized harmonic potential must be finite and nonnegative.');
end
end

function problem = makeProblem(grid, V, template, potentialCaseLabel)
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
problem.potential_label = potentialCaseLabel;
problem.entropy.enabled = false;
problem.entropy.eta = 0;
end

function rates = rateSeries(errors, floorValue)
errors = errors(:).';
rates = NaN(size(errors));
for index = 2:numel(errors)
    if errors(index - 1) > floorValue && errors(index) > floorValue
        rates(index) = log(errors(index - 1) / errors(index)) / log(2);
    end
end
end

function accelerated = detectAcceleration(records)
rates = [records.rate_L2];
finiteRates = rates(isfinite(rates));
if numel(finiteRates) < 3
    accelerated = false;
else
    accelerated = max(finiteRates) >= 4 ...
        && max(finiteRates(end - 1:end)) ...
        >= median(finiteRates(1:min(2, end))) + 1.5;
end
end

function printMainTable(records, epsilon, sigma, L, R0, R1, Nref, floorValue)
fprintf('\n============================================================\n');
fprintf(['Fourier mesh refinement\n' ...
    'harmonic trap with C-infinity boundary periodicization\n']);
fprintf('============================================================\n');
fprintf('epsilon = %.3e\n', epsilon);
fprintf('sigma   = %.3e\n', sigma);
fprintf('L       = %.6g\n', L);
fprintf('R0      = %.6g\n', R0);
fprintf('R1      = %.6g\n', R1);
fprintf('Nref    = %d\n\n', Nref);
fprintf(' N       L2_total   rate      Linf_res   rate      dE          rate      PG_final\n');
for index = 1:numel(records)
    r = records(index);
    fprintf('%4d   %.3e  %7s   %.3e  %7s   %.3e  %7s   %.3e\n', ...
        r.N, r.total_L2, rateLabel(r.rate_L2, r.total_L2, floorValue), ...
        r.resolved_Linf, rateLabel(r.rate_Linf, r.resolved_Linf, floorValue), ...
        r.target_energy_error, rateLabel(r.rate_energy, ...
        r.target_energy_error, floorValue), r.final_pg);
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

function label = passFail(value)
if value
    label = 'YES';
else
    label = 'NO';
end
end

function writeLatexTable(records, filename, floorValue)
file = fopen(filename, 'w');
if file < 0
    error('Unable to create %s.', filename);
end
cleanup = onCleanup(@() fclose(file));
fprintf(file, '\\begin{tabular}{rcccccc}\n');
fprintf(file, '\\hline\n');
fprintf(file, ['$N$ & $L^2$ error & rate & $L^\\infty$ error & ' ...
    'rate & energy error & rate \\\\\n']);
fprintf(file, '\\hline\n');
for index = 1:numel(records)
    r = records(index);
    fprintf(file, '%d & %.6e & %s & %.6e & %s & %.6e & %s \\\\\n', ...
        r.N, r.total_L2, latexRate(r.rate_L2, r.total_L2, floorValue), ...
        r.resolved_Linf, latexRate(r.rate_Linf, r.resolved_Linf, floorValue), ...
        r.target_energy_error, latexRate(r.rate_energy, ...
        r.target_energy_error, floorValue));
end
fprintf(file, '\\hline\n');
fprintf(file, '\\end{tabular}\n');
clear cleanup;
end

function label = latexRate(rate, errorValue, floorValue)
if errorValue <= floorValue
    label = '--';
elseif isfinite(rate)
    label = sprintf('%.2f', rate);
else
    label = '--';
end
end

function makeFigures(solutions, records, Vperiodic, Vharmonic, ...
    referenceGrid, R0, R1, folder, floorValue)
Nvalues = [records.N];
f1 = figure('Visible', 'off', 'Name', 'Spectral harmonic state convergence');
loglog(Nvalues, max([records.total_L2], floorValue), 'o-', ...
    Nvalues, max([records.resolved_Linf], floorValue), 's-');
grid on; xlabel('N'); ylabel('state error');
legend('L2 total', 'Linf resolved', 'Location', 'best');
exportgraphics(f1, fullfile(folder, '01_state_convergence.png'), ...
    'Resolution', 180);
close(f1);

f2 = figure('Visible', 'off', 'Name', 'Spectral harmonic energy convergence');
loglog(Nvalues, max([records.target_energy_error], floorValue), 'o-');
grid on; xlabel('N'); ylabel('|E_N-E_{ref}|');
exportgraphics(f2, fullfile(folder, '02_energy_convergence.png'), ...
    'Resolution', 180);
close(f2);

f3 = figure('Visible', 'off', 'Name', 'Spectral harmonic Fourier decay');
hold on;
representative = [64, 128, 256, 512, referenceGrid.N];
allN = cellfun(@(state) state.grid.N, solutions);
for N = representative
    index = find(allN == N, 1);
    [coefficients, modes] = ...
        src.discretization.ps.FourierCoefficients(solutions{index}.rho);
    [absoluteModes, order] = sort(abs(modes));
    semilogy(absoluteModes, max(abs(coefficients(order)), realmin), ...
        '.', 'DisplayName', sprintf('N=%d', N));
end
grid on; xlabel('|k|'); ylabel('|rho hat_k|');
legend('Location', 'best');
exportgraphics(f3, fullfile(folder, '03_fourier_decay.png'), ...
    'Resolution', 180);
close(f3);

f4 = figure('Visible', 'off', 'Name', 'Harmonic boundary periodicization');
mask = referenceGrid.x >= 0;
plot(referenceGrid.x(mask), Vharmonic(mask), '-', ...
    referenceGrid.x(mask), Vperiodic(mask), '--', 'LineWidth', 1.2);
hold on; xline(R0, ':', 'R0'); xline(R1, ':', 'R1');
grid on; xlabel('x'); ylabel('V(x)');
legend('x^2/2', 'C-infinity periodicized harmonic', ...
    'Location', 'best');
exportgraphics(f4, fullfile(folder, ...
    '04_harmonic_periodicization.png'), 'Resolution', 180);
close(f4);
end
