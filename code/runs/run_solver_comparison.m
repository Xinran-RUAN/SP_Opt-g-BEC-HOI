%RUN_SOLVER_COMPARISON Compare first-order solvers with one shared PDAS polish.
clearvars; clc;

% -------------------------- editable settings --------------------------
regularization_name = 'shift_smooth';
epsilon = 1e-3;
N = 256;
L = 32;
beta = 10;
delta = 10;

solver_names = {'ista', 'spg', 'fista_cd'};
projection_name = 'simplex';
first_stage_pg_tol = 1e-8;
polish_pg_tol = 1e-12;
max_iter = 200000;

show_plot = true;
save_result = true;
overwrite_existing = false;
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
baseConfig = experiments.DefaultConfig();
baseConfig.parameters.beta = beta;
baseConfig.parameters.delta = delta;
baseConfig.parameters.L = L;
baseConfig.parameters.N = N;
baseConfig.parameters.epsilon = epsilon;
baseConfig.regularization.name = regularization_name;
baseConfig.regularization.epsilon = epsilon;
baseConfig.regularization.transition_width = epsilon;
baseConfig.solver.projection_name = projection_name;
baseConfig.solver.pg_tol = first_stage_pg_tol;
baseConfig.solver.polish_mode = 'always';
baseConfig.solver.polish.pg_tol = polish_pg_tol;
baseConfig.solver.max_iter = max_iter;
baseConfig.solver.display = false;

commonGrid = model.SetupGrid1D(baseConfig.parameters);
rho0 = model.InitialDensity(commonGrid, ...
    baseConfig.parameters.mass, projection_name);
numberOfSolvers = numel(solver_names);
records = repmat(emptyRecord(), numberOfSolvers, 1);
rho = cell(numberOfSolvers, 1);
energy = nan(numberOfSolvers, 1);
history = cell(numberOfSolvers, 1);

for j = 1:numberOfSolvers
    config = baseConfig;
    config.solver.name = solver_names{j};
    fprintf('\nsolver comparison: %s\n', upper(strrep(solver_names{j}, '_', '-')));
    result = experiments.SolveGroundState(config, rho0);
    records(j) = recordFromResult(solver_names{j}, result);
    rho{j} = result.rho;
    energy(j) = result.energy;
    history{j} = result.history;
    fprintf(['  main iter=%d time=%.3fs E=%.12e PG=%.3e; ' ...
        'polish iter=%d time=%.3fs; Efinal=%.12e PG=%.3e KKT=%.3e\n'], ...
        records(j).first_stage_iterations, records(j).first_stage_time, ...
        records(j).first_stage_energy, records(j).first_stage_pg_residual, ...
        records(j).polish_iterations, records(j).polish_time, ...
        records(j).final_energy, records(j).final_pg_residual, ...
        records(j).final_kkt_residual);
end

parameters = baseConfig.parameters;
regularization.name = regularization_name;
regularization.epsilon = epsilon;
regularization.transition_width = epsilon;
grid = commonGrid;
solver = baseConfig.solver;
diagnostics.records = records;
diagnostics.solver_names = solver_names;
archive.parameters = parameters;
archive.regularization = regularization;
archive.grid = grid;
archive.solver = solver;
archive.rho = rho;
archive.energy = energy;
archive.diagnostics = diagnostics;
archive.history = history;

if save_result
    baseName = sprintf('solver_comparison_N%d_eps_%s.mat', ...
        N, strrep(sprintf('%.0e', epsilon), '-', 'm'));
    result_file = experiments.SaveResult( ...
        archive, 'diagnostics', baseName, overwrite_existing);
    if ~isempty(result_file)
        fprintf('saved: %s\n', result_file);
    end
end
if show_plot && save_result && exist('result_file', 'var') && ~isempty(result_file)
    run(fullfile(root, 'post', 'post_solver_comparison.m'));
end

function record = emptyRecord()
record.solver_name = '';
record.first_stage_iterations = NaN;
record.first_stage_time = NaN;
record.first_stage_energy = NaN;
record.first_stage_pg_residual = NaN;
record.polish_iterations = NaN;
record.polish_time = NaN;
record.final_energy = NaN;
record.final_pg_residual = NaN;
record.final_kkt_residual = NaN;
end

function record = recordFromResult(solverName, result)
record = emptyRecord();
record.solver_name = solverName;
record.first_stage_iterations = result.diagnostics.main_iterations;
record.first_stage_time = result.diagnostics.main_elapsed_time;
record.first_stage_energy = result.diagnostics.energy_before_polish;
record.first_stage_pg_residual = result.diagnostics.pg_before_polish;
record.polish_iterations = result.diagnostics.polish_iterations;
record.polish_time = result.diagnostics.polish_elapsed_time;
record.final_energy = result.energy;
record.final_pg_residual = result.diagnostics.pg_residual;
record.final_kkt_residual = result.diagnostics.kkt_residual;
end
