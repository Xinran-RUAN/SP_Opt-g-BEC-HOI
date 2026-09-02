%RUN_ENTROPY_MESH_REFINEMENT Fixed-eta spectral-space refinement study.
clearvars; clc;

% -------------------------- editable settings --------------------------
eta = 1e-1;
epsilon = 1e-3;
N_list = [32, 64, 128, 256, 512];
N_ref = 1024;
energy_diag_factor = 2;
L = 32;
beta = 10;
delta = 10;
pg_tol = 1e-8;
polish_kkt_tol = 1e-10;
reference_pg_tol = 1e-10;
reference_tail_tol = 1e-12;
max_iter = 200000;
show_plot = true;
save_result = true;
overwrite_existing = true;
% -----------------------------------------------------------------------

if eta <= 0 || any(mod(N_list, 2) ~= 0) || mod(N_ref, 2) ~= 0 ...
        || N_ref <= max(N_list) || any(mod(N_ref, N_list) ~= 0)
    error(['Use eta>0 and strictly nested even meshes with ' ...
        'N_ref>max(N_list).']);
end
N_diag = energy_diag_factor * N_ref;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
config = experiments.DefaultConfig();
config.parameters.beta = beta;
config.parameters.delta = delta;
config.parameters.L = L;
config.parameters.epsilon = epsilon;
config.regularization.name = 'shift_smooth';
config.regularization.epsilon = epsilon;
config.regularization.transition_width = epsilon;
config.entropy.enabled = true;
config.entropy.eta = eta;
config.solver.name = 'fista_cd';
config.solver.pg_tol = pg_tol;
config.solver.max_iter = max_iter;
config.solver.polish_mode = 'always';
config.solver.polish.pg_tol = polish_kkt_tol;
config.solver.display = false;

allN = [N_list, N_ref];
solutions = cell(numel(allN), 1);
warmStart = [];
transferProblem = [];
for j = 1:numel(allN)
    config.parameters.N = allN(j);
    fineGrid = model.SetupGrid1D(config.parameters);
    if ~isempty(warmStart)
        transferProblem.grid = fineGrid;
        transferProblem.plan = src.discretization.ps.Plan1D(fineGrid);
        warmStart = experiments.TransferState1D( ...
            warmStart, fineGrid, config.parameters.mass, ...
            transferProblem, config.solver);
    end
    solutions{j} = experiments.SolveGroundState(config, warmStart);
    warmStart = solutions{j}.rho;
    transferProblem = solutions{j}.problem;
end

finestState = solutions{end};
reference_ok = finestState.diagnostics.pg_residual <= reference_pg_tol ...
    && finestState.diagnostics.fourier_tail.tail_ratio_quarter ...
    <= reference_tail_tol;
if ~reference_ok
    warning(['N_ref state is not sufficiently resolved; reported errors ' ...
        'are finest-grid differences only.']);
    comparisonLabel = 'finest-grid comparison state';
else
    comparisonLabel = 'reference state';
end
commonFinest = src.diagnostics.CommonGridEnergy( ...
    finestState.rho, finestState.grid, finestState.problem, N_diag);
firstRecord = src.diagnostics.MeshRefinementComparison( ...
    solutions{1}, finestState, N_diag, commonFinest);
records = repmat(firstRecord, numel(N_list), 1);
for j = 2:numel(N_list)
    records(j) = src.diagnostics.MeshRefinementComparison( ...
        solutions{j}, finestState, N_diag, commonFinest);
end

fprintf('\nfixed eta=%.3e spectral comparison against %s N=%d\n', ...
    eta, comparisonLabel, N_ref);
fprintf('N      Eaug(native)   L2-res      L2-tail     L2-total    Linf-res    PG/KKT               FFT-tail    regime\n');
for j = 1:numel(records)
    fprintf('%4d   %.3e      %.3e   %.3e   %.3e   %.3e   %.2e/%.2e   %.3e   %s\n', ...
        records(j).N, records(j).native_augmented_energy_error, ...
        records(j).resolved_L2_error, records(j).reference_tail_L2, ...
        records(j).total_spectral_L2_error, ...
        records(j).resolved_Linf_error, records(j).pg_residual, ...
        records(j).kkt_residual, records(j).fft_tail_ratio, ...
        records(j).regime);
end
fprintf('\nN      dEphys-native   dEphys-common   common phys/aug valid   tail-mass    tail-max\n');
for j = 1:numel(records)
    fprintf('%4d   %+13.3e   %13.3e   %d/%d                    %.3e   %.3e\n', ...
        records(j).N, records(j).native_physical_energy_difference, ...
        records(j).common_physical_energy_error, ...
        records(j).common_physical_energy_valid, ...
        records(j).common_augmented_energy_valid, ...
        records(j).tail_mass, records(j).tail_max);
end
if finestState.diagnostics.tail_indicator > 1e-8
    warning('Increase L before interpreting Fourier convergence.');
end

archive.parameters = config.parameters;
archive.parameters.N_list = N_list;
archive.parameters.N_ref = N_ref;
archive.parameters.N_diag = N_diag;
archive.parameters.eta = eta;
archive.regularization = finestState.regularization;
archive.grid = cellfun(@(s) s.grid, solutions, 'UniformOutput', false);
archive.solver = config.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = cellfun(@(s) s.physical_energy, solutions);
archive.physical_energy = archive.energy;
archive.entropy_value = cellfun(@(s) s.entropy_value, solutions);
archive.augmented_energy = cellfun(@(s) s.augmented_energy, solutions);
archive.diagnostics.records = records;
archive.diagnostics.finest = finestState.diagnostics;
archive.diagnostics.reference_ok = reference_ok;
archive.diagnostics.comparison_label = comparisonLabel;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'entropy_mesh', ...
        sprintf('entropy_mesh_eta_%s_Nref%d.mat', ...
        strrep(sprintf('%.0e', eta), '-', 'm'), N_ref), overwrite_existing);
end
if show_plot
    makePlots(solutions, records, N_list, root, save_result, eta);
end

function makePlots(solutions, records, Nvalues, root, saveFigures, eta)
folder = fullfile(root, 'results', 'entropy_mesh');
if saveFigures && ~isfolder(folder)
    mkdir(folder);
end
prefix = sprintf('eta_%s_', strrep(sprintf('%.0e', eta), '-', 'm'));

f1 = figure('Name', 'Spectral state errors');
loglog(Nvalues, [records.resolved_L2_error], 'o-', ...
    Nvalues, [records.reference_tail_L2], 's-', ...
    Nvalues, [records.total_spectral_L2_error], '^-');
xlabel('N'); ylabel('L2 error'); grid on;
legend('resolved', 'reference tail', 'total', 'Location', 'best');
saveFigure(f1, folder, [prefix '01_spectral_state.png'], saveFigures);

f2 = figure('Name', 'Energy diagnostics');
loglog(Nvalues, [records.native_augmented_energy_error], 'o-', ...
    Nvalues, [records.common_physical_energy_error], 's-');
xlabel('N'); ylabel('energy diagnostic'); grid on;
legend('native augmented', 'common physical', 'Location', 'best');
saveFigure(f2, folder, [prefix '02_energy.png'], saveFigures);

f3 = figure('Name', 'Fourier coefficient decay'); hold on;
representative = unique([1, ceil(numel(Nvalues)/2), numel(Nvalues)]);
for q = representative
    [coefficients, modes] = ...
        src.discretization.ps.FourierCoefficients(solutions{q}.rho);
    semilogy(abs(modes), abs(coefficients), '.', ...
        'DisplayName', sprintf('N=%d', Nvalues(q)));
end
set(gca, 'YScale', 'log'); xlabel('|k|'); ylabel('abs(rho hat)');
grid on; legend('Location', 'best');
saveFigure(f3, folder, [prefix '03_fourier.png'], saveFigures);

f4 = figure('Name', 'Optimization versus spatial error');
loglog(Nvalues, [records.total_spectral_L2_error], 'o-', ...
    Nvalues, [records.pg_residual], 's-', ...
    Nvalues, [records.kkt_residual], '^-');
xlabel('N'); ylabel('residual or state error'); grid on;
legend('total spectral L2', 'PG', 'KKT', 'Location', 'best');
saveFigure(f4, folder, [prefix '04_optimization_floor.png'], saveFigures);
end

function saveFigure(handle, folder, name, enabled)
if enabled
    exportgraphics(handle, fullfile(folder, name), 'Resolution', 180);
end
end
