%RUN_DIAG_PERIODIC_POTENTIAL_MESH Optional non-production mesh experiment.
%
% Run manually only after run_spatial_operator_consistency identifies the
% non-smooth periodic extension of harmonic V as the obstruction. This file
% changes only the diagnostic problem's V; it does not modify production
% potential construction or solver settings.
clearvars; clc;

epsilon = 1e-3;
sigma = 1e-12;
beta = 10;
delta = 10;
mass = 1;
L = 8;
N_list = [32, 64, 128, 256, 512];
N_ref = 1024;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
base = experiments.DefaultConfig();
solver = base.solver;
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
solver.display = true;
solver.display_every = 200;

allN = [N_list, N_ref];
solutions = cell(numel(allN), 1);
rho0 = [];
for j = 1:numel(allN)
    grid = makeGrid(L, allN(j));
    problem = makeProblem(grid, epsilon, sigma, beta, delta, mass, L);
    if isempty(rho0)
        rho0 = exp(-grid.x .^ 2) / sqrt(pi);
    else
        rho0 = src.discretization.ps.Prolong(rho0, grid.N);
    end
    rho0 = src.constraints.ProjectPositiveConservative( ...
        rho0, mass, grid.h, 1e-14);
    result = src.SolveGroundState1D(problem, rho0, solver);
    result.grid = grid;
    result.fourier = src.diagnostics.FourierTailDiagnostics(result.rho);
    solutions{j} = result;
    rho0 = result.rho;
end

reference = solutions{end};
records = repmat(struct(), numel(N_list), 1);
for j = 1:numel(N_list)
    comparison = src.diagnostics.SpectralStateComparison( ...
        solutions{j}.rho, solutions{j}.grid, ...
        reference.rho, reference.grid);
    records(j).N = N_list(j);
    records(j).resolved_L2 = comparison.resolved_L2_error;
    records(j).reference_tail_L2 = comparison.reference_tail_L2;
    records(j).total_L2 = comparison.total_spectral_L2_error;
    records(j).target_energy_error = abs( ...
        solutions{j}.target_energy - reference.target_energy);
    records(j).final_pg = solutions{j}.diagnostics.final_pg_residual;
end

fprintf('\nDIAGNOSTIC smooth-periodic-potential minimizer mesh\n');
fprintf(' N       L2_res       L2_tail      L2_total     dE_target     finalPG\n');
for j = 1:numel(records)
    r = records(j);
    fprintf('%4d  %.3e  %.3e  %.3e  %.3e  %.3e\n', ...
        r.N, r.resolved_L2, r.reference_tail_L2, r.total_L2, ...
        r.target_energy_error, r.final_pg);
end

archive.description = ['Diagnostic smooth periodic V only; ' ...
    'not a production-model result.'];
archive.N_list = N_list;
archive.N_ref = N_ref;
archive.records = records;
archive.rho = cellfun(@(s) s.rho, solutions, 'UniformOutput', false);
archive.grid = cellfun(@(s) s.grid, solutions, 'UniformOutput', false);
archive.energy = cellfun(@(s) s.target_energy, solutions);
archive.diagnostics = cellfun(@(s) s.diagnostics, solutions, ...
    'UniformOutput', false);
outputFolder = fullfile(root, 'results', 'diagnostics');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
save(fullfile(outputFolder, 'DIAG_periodic_potential_mesh.mat'), ...
    'archive', '-v7.3');

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end

function problem = makeProblem(grid, epsilon, sigma, beta, delta, mass, L)
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = (L ^ 2 / pi ^ 2) * (1 - cos(pi * grid.x / L));
problem.beta = beta;
problem.delta = delta;
problem.mass = mass;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(rho) rho + epsilon;
problem.fisher_regularization.ds_epsilon = @(rho) ones(size(rho));
problem.fisher_regularization.d2s_epsilon = @(rho) zeros(size(rho));
problem.fisher_regularization.label = 's_epsilon(rho) = rho + epsilon';
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(rho) ...
    rho ./ hypot(rho, sigma);
problem.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
problem.potential_regularization.label = ...
    'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
problem.potential_regularization.prox_type = 'generic_convex';
problem.entropy.enabled = false;
problem.entropy.eta = 0;
end
