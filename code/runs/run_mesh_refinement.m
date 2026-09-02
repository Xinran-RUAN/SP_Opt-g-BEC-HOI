%RUN_MESH_REFINEMENT eta=0 spectral-space refinement diagnostics.
clearvars; clc;

% -------------------------- editable settings --------------------------
regularization_names = {
    'shift_smooth'
    'piecewise_c1'
    'piecewise_c2'
    'piecewise_c3'
};
epsilon = 1e-3;
transition_width = epsilon;
N_list = [32, 64, 128, 256, 512];
N_ref = 1024;
energy_diag_factor = 2;
L = 32;
beta = 10;
delta = 10;
solver_name = 'spg';
projection_name = 'simplex';
pg_tol = 1e-7;
max_iter = 200000;
polish_mode = 'if_needed';
polish_pg_tol = 1e-12;
reference_pg_tol = 1e-10;
reference_tail_tol = 1e-12;
show_plot = true;
save_result = true;
overwrite_existing = false;
% -----------------------------------------------------------------------

if any(mod(N_list, 2) ~= 0) || mod(N_ref, 2) ~= 0 ...
        || N_ref <= max(N_list) || any(mod(N_ref, N_list) ~= 0)
    error('Use strictly nested even meshes with N_ref>max(N_list).');
end
N_diag = energy_diag_factor * N_ref;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
baseConfig = experiments.DefaultConfig();
baseConfig.parameters.beta = beta;
baseConfig.parameters.delta = delta;
baseConfig.parameters.L = L;
baseConfig.parameters.epsilon = epsilon;
baseConfig.solver.name = solver_name;
baseConfig.solver.projection_name = projection_name;
baseConfig.solver.pg_tol = pg_tol;
baseConfig.solver.max_iter = max_iter;
baseConfig.solver.polish_mode = polish_mode;
baseConfig.solver.polish.pg_tol = polish_pg_tol;
baseConfig.solver.display = false;
baseConfig.entropy.enabled = false;
baseConfig.entropy.eta = 0;

numberOfRegularizations = numel(regularization_names);
numberOfMeshes = numel(N_list);
solutions = cell(numberOfRegularizations, numberOfMeshes + 1);
records = cell(numberOfRegularizations, 1);
reference_ok = false(numberOfRegularizations, 1);
comparison_labels = cell(numberOfRegularizations, 1);
regularization_info = cell(numberOfRegularizations, 1);

for regIndex = 1:numberOfRegularizations
    regName = regularization_names{regIndex};
    config = baseConfig;
    config.regularization.name = regName;
    config.regularization.epsilon = epsilon;
    config.regularization.transition_width = transition_width;
    warmStart = [];
    fprintf('\nmesh refinement: %s\n', regName);
    for meshIndex = 1:numberOfMeshes + 1
        if meshIndex <= numberOfMeshes
            currentN = N_list(meshIndex);
        else
            currentN = N_ref;
        end
        config.parameters.N = currentN;
        fineGrid = model.SetupGrid1D(config.parameters);
        if ~isempty(warmStart)
            warmStart = experiments.TransferState1D( ...
                warmStart, fineGrid, config.parameters.mass);
        end
        solutions{regIndex, meshIndex} = ...
            experiments.SolveGroundState(config, warmStart);
        warmStart = solutions{regIndex, meshIndex}.rho;
    end

    finestState = solutions{regIndex, end};
    regularization_info{regIndex} = finestState.regularization;
    reference_ok(regIndex) = ...
        finestState.diagnostics.pg_residual <= reference_pg_tol ...
        && finestState.diagnostics.fourier_tail.tail_ratio_quarter ...
        <= reference_tail_tol;
    if reference_ok(regIndex)
        comparison_labels{regIndex} = 'reference state';
    else
        comparison_labels{regIndex} = 'finest-grid comparison state';
        warning(['%s N_ref state is not sufficiently resolved; reported ' ...
            'errors are finest-grid differences only.'], regName);
    end
    commonFinest = src.diagnostics.CommonGridEnergy( ...
        finestState.rho, finestState.grid, finestState.problem, N_diag);
    firstRecord = src.diagnostics.MeshRefinementComparison( ...
        solutions{regIndex, 1}, finestState, N_diag, commonFinest);
    regRecords = repmat(firstRecord, numberOfMeshes, 1);
    for meshIndex = 2:numberOfMeshes
        regRecords(meshIndex) = src.diagnostics.MeshRefinementComparison( ...
            solutions{regIndex, meshIndex}, finestState, N_diag, commonFinest);
    end
    records{regIndex} = regRecords;

    fprintf('comparison against %s N=%d\n', ...
        comparison_labels{regIndex}, N_ref);
    fprintf('N      Eaug(native)   L2-res      L2-tail     L2-total    Linf-res    PG/KKT               FFT-tail    regime\n');
    for j = 1:numberOfMeshes
        r = regRecords(j);
        fprintf('%4d   %.3e      %.3e   %.3e   %.3e   %.3e   %.2e/%.2e   %.3e   %s\n', ...
            r.N, r.native_augmented_energy_error, r.resolved_L2_error, ...
            r.reference_tail_L2, r.total_spectral_L2_error, ...
            r.resolved_Linf_error, r.pg_residual, r.kkt_residual, ...
            r.fft_tail_ratio, r.regime);
    end
    fprintf('N      dEphys-native   dEphys-common   common-valid   tail-mass\n');
    for j = 1:numberOfMeshes
        r = regRecords(j);
        fprintf('%4d   %+13.3e   %13.3e   %d             %.3e\n', ...
            r.N, r.native_physical_energy_difference, ...
            r.common_physical_energy_error, ...
            r.common_physical_energy_valid, r.tail_mass);
    end
end

archive.parameters = baseConfig.parameters;
archive.parameters.N_list = N_list;
archive.parameters.N_ref = N_ref;
archive.parameters.N_diag = N_diag;
archive.regularization.names = regularization_names;
archive.regularization.epsilon = epsilon;
archive.regularization.transition_width = transition_width;
archive.regularization.info = regularization_info;
archive.grid = cellfun(@(s) s.grid, solutions, 'UniformOutput', false);
archive.solver = baseConfig.solver;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.energy = cellfun(@(s) s.physical_energy, solutions);
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.diagnostics.reference_ok = reference_ok;
archive.diagnostics.comparison_labels = comparison_labels;
archive.history = cellfun(@(s) s.history, solutions, 'UniformOutput', false);

if save_result
    baseName = sprintf('mesh_refinement_Nref%d_eps_%s.mat', ...
        N_ref, strrep(sprintf('%.0e', epsilon), '-', 'm'));
    result_file = experiments.SaveResult( ...
        archive, 'mesh_refinement', baseName, overwrite_existing);
end
if show_plot
    makePlots(solutions, records, regularization_names, N_list, ...
        root, save_result);
end

function makePlots(solutions, records, names, Nvalues, root, saveFigures)
folder = fullfile(root, 'results', 'mesh_refinement');
if saveFigures && ~isfolder(folder)
    mkdir(folder);
end
numberOfRegularizations = numel(names);
f1 = figure('Name', 'Spectral state convergence');
tiledlayout(numberOfRegularizations, 1);
f2 = figure('Name', 'Energy convergence');
tiledlayout(numberOfRegularizations, 1);
f3 = figure('Name', 'Fourier coefficient decay');
tiledlayout(numberOfRegularizations, 1);
f4 = figure('Name', 'Optimization versus spatial error');
tiledlayout(numberOfRegularizations, 1);
for regIndex = 1:numberOfRegularizations
    r = records{regIndex};
    figure(f1); nexttile;
    loglog(Nvalues, [r.resolved_L2_error], 'o-', ...
        Nvalues, [r.reference_tail_L2], 's-', ...
        Nvalues, [r.total_spectral_L2_error], '^-');
    grid on; title(names{regIndex}, 'Interpreter', 'none');
    if regIndex == 1
        legend('resolved', 'reference tail', 'total', 'Location', 'best');
    end

    figure(f2); nexttile;
    loglog(Nvalues, [r.native_augmented_energy_error], 'o-', ...
        Nvalues, [r.common_physical_energy_error], 's-');
    grid on; title(names{regIndex}, 'Interpreter', 'none');
    if regIndex == 1
        legend('native augmented', 'common physical', 'Location', 'best');
    end

    figure(f3); nexttile; hold on;
    representative = unique([1, ceil(numel(Nvalues)/2), numel(Nvalues)]);
    for q = representative
        [coefficients, modes] = src.discretization.ps.FourierCoefficients( ...
            solutions{regIndex, q}.rho);
        semilogy(abs(modes), abs(coefficients), '.', ...
            'DisplayName', sprintf('N=%d', Nvalues(q)));
    end
    set(gca, 'YScale', 'log'); grid on;
    title(names{regIndex}, 'Interpreter', 'none');
    if regIndex == 1
        legend('Location', 'best');
    end

    figure(f4); nexttile;
    loglog(Nvalues, [r.total_spectral_L2_error], 'o-', ...
        Nvalues, [r.pg_residual], 's-', Nvalues, [r.kkt_residual], '^-');
    grid on; title(names{regIndex}, 'Interpreter', 'none');
    if regIndex == 1
        legend('total L2', 'PG', 'KKT', 'Location', 'best');
    end
end
if saveFigures
    exportgraphics(f1, fullfile(folder, '01_spectral_state.png'), 'Resolution', 180);
    exportgraphics(f2, fullfile(folder, '02_energy.png'), 'Resolution', 180);
    exportgraphics(f3, fullfile(folder, '03_fourier.png'), 'Resolution', 180);
    exportgraphics(f4, fullfile(folder, '04_optimization_floor.png'), 'Resolution', 180);
end
end
