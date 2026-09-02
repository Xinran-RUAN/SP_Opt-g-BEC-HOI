%RUN_POLISH_LINEAR_SOLVER_COMPARISON Same-handoff GMRES versus PCG-Schur.
clearvars; clc;

potential_name = 'sqrt_squared_scale';
epsilon = 1e-3;
L = 8;
N_list = [128, 256];
save_result = true;
overwrite_existing = true;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
baseConfig = experiments.DefaultConfig();
baseConfig.parameters.beta = 10;
baseConfig.parameters.delta = 10;
baseConfig.parameters.mass = 1;
baseConfig.parameters.L = L;
baseConfig.parameters.epsilon = epsilon;
baseConfig.regularization.name = 'shift_smooth';
baseConfig.regularization.epsilon = epsilon;
baseConfig.regularization.transition_width = epsilon;
baseConfig.potential_regularization.name = potential_name;
baseConfig.entropy.enabled = false;
baseConfig.entropy.eta = 0;
baseConfig.solver.name = 'fista_cd';
baseConfig.solver.splitting = 'potential_prox';
baseConfig.solver.projection_name = 'semismooth';
baseConfig.solver.projection_tol = 1e-14;
baseConfig.solver.residual_check_interval = 10;
baseConfig.solver.max_iter = 20000;
baseConfig.solver.switch.enabled = true;
baseConfig.solver.switch.energy_window = 50;
baseConfig.solver.switch.energy_tol = 1e-12;
baseConfig.solver.switch.consecutive_windows = 2;
baseConfig.solver.switch.min_iter = 200;
baseConfig.solver.switch.pg_entry_tol = 1e-5;
baseConfig.solver.switch.max_main_iter = 20000;
baseConfig.solver.polish_mode = 'none';
baseConfig.solver.display = false;

records = repmat(emptyRecord(), 2 * numel(N_list), 1);
handoffResults = cell(numel(N_list), 1);
polishResults = cell(numel(N_list), 2);
differences = repmat(struct(), numel(N_list), 1);
row = 0;
for j = 1:numel(N_list)
    config = baseConfig;
    config.parameters.N = N_list(j);
    grid = model.SetupGrid1D(config.parameters);
    rho0 = exp(-grid.x .^ 2) / sqrt(pi);
    rho0 = rho0 * (config.parameters.mass ...
        / src.constraints.Mass(rho0, grid.h));
    handoffResults{j} = experiments.SolveGroundState(config, rho0);
    rhoHandoff = handoffResults{j}.rho;
    problem = handoffResults{j}.problem;
    options = handoffResults{j}.solver.polish;
    options.entry_pg_tol = 1e-4;
    options.pg_tol = 1e-12;
    options.max_iter = 20;
    options.pdas_max_iter = 20;
    options.display = false;

    oldOptions = options;
    oldOptions.linear_solver = 'pdas_gmres';
    polishResults{j, 1} = src.solvers.PolishKKT( ...
        rhoHandoff, problem, oldOptions);
    newOptions = options;
    newOptions.linear_solver = 'interior_pcg_schur';
    newOptions.preconditioner = 'fd_variable';
    newOptions.allow_pdas_fallback = true;
    polishResults{j, 2} = src.solvers.PolishKKT( ...
        rhoHandoff, problem, newOptions);

    for method = 1:2
        row = row + 1;
        records(row) = makeRecord( ...
            N_list(j), polishResults{j, method});
    end
    old = polishResults{j, 1};
    new = polishResults{j, 2};
    differences(j).N = N_list(j);
    differences(j).state_L2 = sqrt(grid.h ...
        * sum((old.rho - new.rho) .^ 2));
    differences(j).target_energy = abs(old.energy - new.energy);
end

fprintf(['\nN    method                 outer  Krylov(z/w or total)  ' ...
    'maxKrylov  time       finalPG      energy             fallback\n']);
for j = 1:numel(records)
    r = records(j);
    if strcmp(r.method, 'PDAS-GMRES')
        krylov = sprintf('%d', r.total_gmres);
    else
        krylov = sprintf('%d/%d', r.total_pcg_z, r.total_pcg_w);
    end
    fprintf('%3d  %-21s %5d  %-20s %9d  %7.3f  %.3e  %.15e  %d\n', ...
        r.N, r.method, r.outer_iterations, krylov, r.max_krylov, ...
        r.elapsed_time, r.final_pg, r.target_energy, r.used_fallback);
end
for j = 1:numel(differences)
    fprintf('N=%d old/new: state L2 %.3e, target energy %.3e\n', ...
        differences(j).N, differences(j).state_L2, ...
        differences(j).target_energy);
end

archive.parameters = baseConfig.parameters;
archive.parameters.N_list = N_list;
archive.regularization.name = 'shift_smooth';
archive.potential_regularization.name = potential_name;
archive.grid = cellfun(@(r) r.grid, handoffResults, ...
    'UniformOutput', false);
archive.solver = baseConfig.solver;
archive.rho = cellfun(@(r) r.rho, polishResults, 'UniformOutput', false);
archive.energy = cellfun(@(r) r.energy, polishResults);
archive.target_energy = archive.energy;
archive.baseline_energy = [];
archive.physical_energy = archive.energy;
archive.entropy_value = zeros(size(archive.energy));
archive.augmented_energy = archive.energy;
archive.diagnostics.records = records;
archive.diagnostics.differences = differences;
archive.history = cellfun(@(r) r.history, polishResults, ...
    'UniformOutput', false);
if save_result
    experiments.SaveResult(archive, 'diagnostics', ...
        'polish_linear_solver_N128_256.mat', overwrite_existing);
end

function record = makeRecord(N, result)
record = emptyRecord();
record.N = N;
if strcmp(result.polish_solver, 'pdas_gmres')
    record.method = 'PDAS-GMRES';
else
    record.method = 'Interior-PCG-Schur';
end
record.outer_iterations = result.iterations;
record.total_gmres = result.total_gmres_iterations;
record.total_pcg_z = result.total_pcg_z_iterations;
record.total_pcg_w = result.total_pcg_w_iterations;
record.max_krylov = max([result.max_gmres_iterations, ...
    result.max_pcg_iterations]);
record.elapsed_time = result.elapsed_time;
record.final_pg = result.pg_residual;
record.target_energy = result.energy;
record.used_fallback = result.switched_to_pdas;
record.status = result.status;
end

function record = emptyRecord()
record.N = NaN;
record.method = '';
record.outer_iterations = NaN;
record.total_gmres = 0;
record.total_pcg_z = 0;
record.total_pcg_w = 0;
record.max_krylov = 0;
record.elapsed_time = NaN;
record.final_pg = NaN;
record.target_energy = NaN;
record.used_fallback = false;
record.status = '';
end
