%RUN_POTENTIAL_MESH_REFINEMENT Fixed-variant spectral mesh diagnostics.
clearvars; clc;
% -----------------------------------------------------------------------
beta = 10;
delta = 10;
mass = 1;
L = 32;
N_list = [32, 64, 128, 256, 512, 1024, 2048];
N_ref = 8192;
main_pg_switch_tol = 1e-6;
final_pg_tol = 1e-12;
reference_pg_tol = 1e-12;
reference_tail_tol = 1e-12;
reference_box_tail_tol = 1e-10; % diagnostic only; not in reference_ok
reference_edge_density_tol = 1e-10; % diagnostic only
reference_edge_derivative_tol = 1e-10; % diagnostic only
max_iter = 200000;
projection_name = 'semismooth';
show_plot = true;
save_result = true;
overwrite_existing = true;
% ======================== regularization ===============================
epsilon = 1e-2;
% Fisher: s_epsilon(rho)
s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho) = rho + epsilon';
% Potential: p_sigma(rho). Set sigma directly.
sigma = 1e-2;
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
p_sigma_label = 'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
potential_name = 'inline_p_sigma'; % metadata only
% =======================================================================

% ===================== trapping potential V(x) =========================
trapping_potential_choice = 'harmonic_cinf_periodic';
% Choices:
%   'harmonic'               : V(x)=x^2/2 on the entire box.
%   'harmonic_cinf_periodic' : unchanged harmonic core with a boundary-only
%                              C-infinity periodic continuation.
R0_fraction = 0.75;
R1_fraction = 0.90;
negligible_density_tol = 1e-10; % diagnostic, not a solver criterion
negligible_mass_tol = 1e-10;    % diagnostic, not a solver criterion

R0 = R0_fraction * L;
R1 = R1_fraction * L;
switch lower(trapping_potential_choice)
    case 'harmonic'
        V = @(x) 0.5 * x .^ 2;
        V_label = 'V(x) = x^2/2';
        boundaryPeriodicization = false;
    case 'harmonic_cinf_periodic'
        V = @(x) model.HarmonicCInfPeriodicPotential(x, L, R0, R1);
        V_label = sprintf([ ...
            'V(x)=x^2/2 for |x|<=%.6g; C-infinity flat continuation ' ...
            'to L^2/2 for |x|>=%.6g'], R0, R1);
        boundaryPeriodicization = true;
    otherwise
        error('Unknown trapping_potential_choice "%s".', ...
            trapping_potential_choice);
end
% =======================================================================

% -----------------------------------------------------------------------

if any(mod(N_list, 2) ~= 0) || mod(N_ref, 2) ~= 0 ...
        || N_ref <= max(N_list) || any(mod(N_ref, N_list) ~= 0)
    error('Use strictly nested even meshes with N_ref>max(N_list).');
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
fprintf('Trapping potential: %s\n', V_label);
config = experiments.DefaultConfig();
config.trapping_potential.V = V;
config.trapping_potential.label = V_label;
config.trapping_potential.mode = trapping_potential_choice;
config.trapping_potential.boundary_periodicized = boundaryPeriodicization;
config.trapping_potential.modification_start = R0;
config.trapping_potential.transition_end = R1;
config.trapping_potential.reference_V = @(x) 0.5 * x .^ 2;
config.trapping_potential.negligible_density_tol = negligible_density_tol;
config.trapping_potential.negligible_mass_tol = negligible_mass_tol;
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.mass = mass;
config.parameters.L = L;
config.parameters.epsilon = epsilon;
config.fisher_regularization.epsilon = epsilon;
config.fisher_regularization.s_epsilon = s_epsilon;
config.fisher_regularization.ds_epsilon = ds_epsilon;
config.fisher_regularization.d2s_epsilon = d2s_epsilon;
config.fisher_regularization.label = fisher_label;
config.potential_regularization.name = potential_name;
config.potential_regularization.sigma = sigma;
config.potential_regularization.p_sigma = p_sigma;
config.potential_regularization.dp_sigma = dp_sigma;
config.potential_regularization.d2p_sigma = d2p_sigma;
config.potential_regularization.label = p_sigma_label;
config.potential_regularization.prox_type = 'generic_convex';
config.entropy.enabled = false;
config.entropy.eta = 0;
config.solver.name = 'fista_cd';
config.solver.splitting = 'potential_prox';
config.solver.projection_name = projection_name;
config.solver.projection_tol = 1e-13;
config.solver.pg_tol = main_pg_switch_tol;
config.solver.final_pg_tol = final_pg_tol;
config.solver.certification_tol = final_pg_tol;
config.solver.residual_check_interval = 10;
% The scalar mass residual accumulates nodal root roundoff across the
% domain.  On coarse, large-box grids (for example L=32, N=32), 1e-14 is
% below that numerical floor.  Use the formally certified mass threshold
% while retaining a tight nodal KKT solve.
config.solver.potential_prox.mass_tol = 1e-12;
config.solver.potential_prox.inner_tol = 1e-14;
config.solver.potential_prox.inner_max_iter = 60;
config.solver.max_iter = max_iter;
config.solver.switch.enabled = true;
config.solver.switch.energy_window = 50;
config.solver.switch.energy_tol = 1e-12;
config.solver.switch.consecutive_windows = 2;
config.solver.switch.min_iter = 200;
config.solver.switch.pg_entry_tol = main_pg_switch_tol;
config.solver.switch.max_main_iter = 20000;
config.solver.polish_mode = 'if_needed';
config.solver.polish.pg_tol = final_pg_tol;
config.solver.polish.max_iter = 20;
config.solver.polish.linear_solver = 'interior_pcg_schur';
config.solver.polish.preconditioner = 'fd_variable';
config.solver.polish.allow_pdas_fallback = true;
config.solver.display = true;
config.solver.display_every = 200;

allN = [N_list, N_ref];
solutions = cell(numel(allN), 1);
warmStart = [];
transferProblem = [];
for j = 1:numel(allN)
    config.parameters.N = allN(j);
    fineGrid = model.SetupGrid1D(config.parameters);
    if isempty(warmStart)
        warmStart = exp(-fineGrid.x .^ 2) / sqrt(pi);
        warmStart = warmStart * (mass ...
            / src.constraints.Mass(warmStart, fineGrid.h));
    else
        transferProblem.grid = fineGrid;
        transferProblem.plan = src.discretization.ps.Plan1D(fineGrid);
        transferSolver = config.solver;
        transferSolver.splitting = 'legacy_full_gradient';
        warmStart = experiments.TransferState1D( ...
            warmStart, fineGrid, mass, transferProblem, transferSolver);
    end
    solutions{j} = experiments.SolveGroundState(config, warmStart);
    warmStart = solutions{j}.rho;
    transferProblem = solutions{j}.problem;
end

finestState = solutions{end};
optimizationPass = finestState.diagnostics.final_pg_residual ...
    <= reference_pg_tol;
fourierPass = finestState.diagnostics.fourier_tail.tail_ratio_quarter ...
    <= reference_tail_tol;
boxPass = max(finestState.diagnostics.far_field.tail_mass, ...
    finestState.diagnostics.far_field.tail_max) <= reference_box_tail_tol;
edgeDensityPass = finestState.diagnostics.far_field.edge_density ...
    <= reference_edge_density_tol;
edgeDerivativePass = finestState.diagnostics.far_field.edge_abs_drho ...
    <= reference_edge_derivative_tol;
reference_ok = optimizationPass && fourierPass;
fprintf('\nReference diagnostics\n');
fprintf('  optimization PG : %s (%.3e)\n', ...
    passFail(optimizationPass), ...
    finestState.diagnostics.final_pg_residual);
fprintf('  Fourier tail    : %s (%.3e)\n', ...
    passFail(fourierPass), ...
    finestState.diagnostics.fourier_tail.tail_ratio_quarter);
fprintf('  box tail        : %s (mass %.3e, max %.3e)\n', ...
    passFail(boxPass), finestState.diagnostics.far_field.tail_mass, ...
    finestState.diagnostics.far_field.tail_max);
fprintf('  edge density    : %s (%.3e) [diagnostic]\n', ...
    passFail(edgeDensityPass), ...
    finestState.diagnostics.far_field.edge_density);
fprintf('  edge derivative : %s (%.3e) [diagnostic]\n', ...
    passFail(edgeDerivativePass), ...
    finestState.diagnostics.far_field.edge_abs_drho);
modification = ...
    finestState.diagnostics.trapping_potential_modification;
if modification.enabled
    fprintf('  harmonic core error : %.3e\n', ...
        modification.core_potential_error);
    fprintf('  modified-region rho : max %.3e, mass %.3e\n', ...
        modification.rho_max, modification.mass);
    if ~modification.density_negligible ...
            || ~modification.mass_negligible
        warning(['The C-infinity continuation starts in a region where ' ...
            'the density is not negligible at the configured diagnostic ' ...
            'tolerances. Increase L or move R0 outward before interpreting.']);
    end
end
if ~reference_ok
    warning(['N_ref is only the finest-grid comparison state, not a ' ...
        'certified spectral reference.']);
    comparisonLabel = 'finest-grid comparison state';
else
    comparisonLabel = 'reference state';
end

firstRecord = makeRecord(solutions{1}, finestState);
records = repmat(firstRecord, numel(N_list), 1);
for j = 2:numel(N_list)
    records(j) = makeRecord(solutions{j}, finestState);
end
fprintf('\n%s (sigma=%.3e) mesh refinement against %s N=%d\n', ...
    potential_name, sigma, comparisonLabel, N_ref);
fprintf(['N      dE_target    dE_base      L2_res      L2_tail     ' ...
    'L2_total    Linf_res    PG_final    mainIter polishIter ' ...
    'FFT_tail    tail_mass  solverOK  regime\n']);
for j = 1:numel(records)
    r = records(j);
    fprintf(['%4d   %.3e   %+.3e   %.3e   %.3e   %.3e   %.3e   ' ...
        '%.2e   %8d %10d %.3e   %.3e   %d         %s\n'], ...
        r.N, r.target_energy_error, r.baseline_energy_difference, ...
        r.resolved_L2_error, r.reference_tail_L2, ...
        r.total_spectral_L2_error, r.resolved_Linf_error, ...
        r.final_pg_residual, r.main_iterations, r.polish_iterations, ...
        r.fft_tail_ratio, r.tail_mass, r.solver_ok, r.regime);
end

archive.parameters = config.parameters;
archive.parameters.N_list = N_list;
archive.parameters.N_ref = N_ref;
archive.trapping_potential = finestState.trapping_potential;
archive.regularization = finestState.regularization;
archive.fisher_regularization = finestState.fisher_regularization;
archive.potential_regularization = finestState.potential_regularization;
archive.grid = cellfun(@(s) s.grid, solutions, 'UniformOutput', false);
archive.solver = config.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = cellfun(@(s) s.target_energy, solutions);
archive.target_energy = archive.energy;
archive.baseline_energy = cellfun(@(s) s.baseline_energy, solutions);
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.diagnostics.reference_ok = reference_ok;
archive.diagnostics.reference_optimization_pass = optimizationPass;
archive.diagnostics.reference_fourier_pass = fourierPass;
archive.diagnostics.reference_box_pass = boxPass;
archive.diagnostics.reference_edge_density_pass = edgeDensityPass;
archive.diagnostics.reference_edge_derivative_pass = edgeDerivativePass;
archive.diagnostics.comparison_label = comparisonLabel;
archive.diagnostics.finest = finestState.diagnostics;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'potential_mesh', ...
        sprintf('potential_mesh_%s_sigma_%s_Nref%d_eps_%s.mat', potential_name, ...
        numberTag(sigma), ...
        N_ref, strrep(sprintf('%.0e', epsilon), '-', 'm')), ...
        overwrite_existing);
end
if show_plot
    makePlots(solutions, records, N_list, root, save_result, ...
        potential_name, sigma);
end

function tag = numberTag(value)
tag = strrep(strrep(sprintf('%.0e', value), '-', 'm'), '+', 'p');
end

function label = passFail(tf)
if tf
    label = 'PASS';
else
    label = 'FAIL';
end
end

function record = makeRecord(resultN, resultRef)
state = src.diagnostics.SpectralStateComparison( ...
    resultN.rho, resultN.grid, resultRef.rho, resultRef.grid);
record.N = resultN.grid.N;
record.target_energy_error = abs( ...
    resultN.target_energy - resultRef.target_energy);
record.baseline_energy_difference = ...
    resultN.baseline_energy - resultRef.baseline_energy;
record.resolved_L2_error = state.resolved_L2_error;
record.resolved_Linf_error = state.resolved_Linf_error;
record.reference_tail_L2 = state.reference_tail_L2;
record.total_spectral_L2_error = state.total_spectral_L2_error;
record.fourier_resolved_coeff_error = state.fourier_resolved_coeff_error;
record.prolonged_L2_error = state.prolonged_L2_error;
record.prolonged_Linf_error = state.prolonged_Linf_error;
record.regime = state.regime;
record.pg_residual = resultN.diagnostics.pg_residual;
record.composite_pg_residual = ...
    resultN.diagnostics.composite_pg_residual;
record.full_pg_residual = resultN.diagnostics.full_pg_residual;
record.main_iterations = resultN.diagnostics.main_iterations;
record.main_pg_residual = resultN.diagnostics.main_pg_residual;
record.polish_iterations = resultN.diagnostics.polish_iterations;
record.final_pg_residual = resultN.diagnostics.final_pg_residual;
record.polish_energy_change = resultN.diagnostics.polish_energy_change;
record.solver_ok = record.final_pg_residual ...
    <= resultN.solver.final_pg_tol;
record.fft_tail_ratio = ...
    resultN.diagnostics.fourier_tail.tail_ratio_quarter;
record.tail_mass = resultN.diagnostics.far_field.tail_mass;
record.tail_max = resultN.diagnostics.far_field.tail_max;
end

function makePlots(solutions, records, Nvalues, root, saveFigures, name, sigma)
folder = fullfile(root, 'results', 'potential_mesh');
if saveFigures && ~isfolder(folder)
    mkdir(folder);
end
f1 = figure('Name', 'Potential spectral state errors');
loglog(Nvalues, [records.resolved_L2_error], 'o-', ...
    Nvalues, [records.reference_tail_L2], 's-', ...
    Nvalues, [records.total_spectral_L2_error], '^-');
grid on; xlabel('N'); ylabel('L2 error');
legend('resolved', 'reference tail', 'total', 'Location', 'best');

f2 = figure('Name', 'Potential energy errors');
loglog(Nvalues, [records.target_energy_error], 'o-', ...
    Nvalues, abs([records.baseline_energy_difference]), 's-');
grid on; xlabel('N'); ylabel('energy diagnostic');
legend('target error', '|baseline difference|', 'Location', 'best');

f3 = figure('Name', 'Potential mesh Fourier decay'); hold on;
representativeN = [64, 128, 256, 512];
for q = 1:numel(representativeN)
    index = find(Nvalues == representativeN(q), 1);
    if ~isempty(index)
        [coefficients, modes] = ...
            src.discretization.ps.FourierCoefficients(solutions{index}.rho);
        semilogy(abs(modes), abs(coefficients), '.', ...
            'DisplayName', sprintf('N=%d', representativeN(q)));
    end
end
set(gca, 'YScale', 'log'); grid on; xlabel('|k|'); ylabel('abs(rho hat)');
legend('Location', 'best');

f4 = figure('Name', 'Potential optimization versus spatial error');
loglog(Nvalues, [records.total_spectral_L2_error], 'o-', ...
    Nvalues, [records.final_pg_residual], 's-');
grid on; xlabel('N'); ylabel('error or residual');
legend('total spectral L2', 'PG', 'Location', 'best');
if saveFigures
    prefix = sprintf('%s_sigma_%s_', name, numberTag(sigma));
    exportgraphics(f1, fullfile(folder, [prefix '01_state.png']), 'Resolution', 180);
    exportgraphics(f2, fullfile(folder, [prefix '02_energy.png']), 'Resolution', 180);
    exportgraphics(f3, fullfile(folder, [prefix '03_fourier.png']), 'Resolution', 180);
    exportgraphics(f4, fullfile(folder, [prefix '04_pg_floor.png']), 'Resolution', 180);
end
end
