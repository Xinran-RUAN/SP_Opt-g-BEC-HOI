%RUN_TWO_DIMENSIONAL_EXAMPLES Two representative 2D density examples.
%
% The production vector optimization path is reused unchanged:
% potential-proximal FISTA-CD -> energy-window handoff -> interior
% Newton-PCG. Only tensor Fourier differentiation and the sparse 2D FD
% preconditioner are dimension-specific.

if ~exist('two_dimensional_smoke_test', 'var')
    two_dimensional_smoke_test = false;
end
if ~exist('do_resolution_check', 'var')
    do_resolution_check = false;
end
if ~exist('reuse_two_dimensional_results', 'var')
    reuse_two_dimensional_results = true;
end
clearvars -except two_dimensional_smoke_test do_resolution_check ...
    reuse_two_dimensional_results;
clc;

% =========================== editable setup ===========================
beta = 10;
delta = 10;
mass = 1;
epsilon = 1e-2;
sigma = 1e-4;
L = 16;
Nx = 512;
Ny = 512;

gamma_x = 1;
gamma_y = 2;
V0 = 5;
lattice_k = pi / 2;
R0_fraction = 0.75;
R1_fraction = 0.90;

% Fisher shift denominator, stated directly in the run file.
r_epsilon = @(rho) rho + epsilon;
dr_epsilon = @(rho) ones(size(rho));
d2r_epsilon = @(rho) zeros(size(rho));

% Stable square-root potential smoothing formulas.
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
% =====================================================================

if two_dimensional_smoke_test
    Nx = 128;
    Ny = 128;
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();
R0 = R0_fraction * L;
R1 = R1_fraction * L;
grid = model.SetupGrid2D(struct('L', L, 'Nx', Nx, 'Ny', Ny));
plan = src.discretization.ps.Plan2D(grid);

[Vx, continuationX] = model.HarmonicCInfPeriodicPotential( ...
    grid.x, L, R0, R1);
[Vy, continuationY] = model.HarmonicCInfPeriodicPotential( ...
    grid.y, L, R0, R1);
V_ah = gamma_x ^ 2 * reshape(Vx, 1, []) ...
    + gamma_y ^ 2 * reshape(Vy, [], 1);
n_period = 2 * lattice_k * L / pi;
if abs(n_period - round(n_period)) > 1e-12
    error('Optical lattice is not periodic on the Fourier box.');
end
V_lattice = V0 * (sin(lattice_k * grid.X) .^ 2 ...
    + sin(lattice_k * grid.Y) .^ 2);
V_optical = V_ah + V_lattice;

coreMask = abs(grid.X) <= R0 & abs(grid.Y) <= R0;
expectedCore = 0.5 * (gamma_x ^ 2 * grid.X .^ 2 ...
    + gamma_y ^ 2 * grid.Y .^ 2);
corePotentialError = max(abs(V_ah(coreMask) - expectedCore(coreMask)));
assert(corePotentialError <= 10 * eps(max(1, max(expectedCore(:)))), ...
    'The anisotropic harmonic potential changed inside its core.');

rho0 = exp(-0.5 * (grid.X .^ 2 + grid.Y .^ 2));
rho0 = rho0(:);
rho0 = mass * rho0 / (grid.h * sum(rho0));
assert(all(rho0 > 0), 'The common Gaussian initial density is not strict.');
assert(abs(grid.h * sum(rho0) - mass) <= 1e-13, ...
    'The common initial density has incorrect mass.');

fisher.r_epsilon = r_epsilon;
fisher.dr_epsilon = dr_epsilon;
fisher.d2r_epsilon = d2r_epsilon;
fisher.epsilon = epsilon;
fisher.label = 'r_epsilon(rho) = rho + epsilon';
potential.p_sigma = p_sigma;
potential.dp_sigma = dp_sigma;
potential.d2p_sigma = d2p_sigma;
potential.sigma = sigma;
potential.label = 'p_sigma(rho)=sqrt(rho^2+sigma^2)-sigma';
potential.name = 'custom_inline';
potential.power = NaN;
potential.prox_type = 'generic_convex';
potential.convexity_tol = 1e-14;
solver = productionSolver();

resultFolder = fullfile(root, 'results', 'two_dimensional_examples');
if two_dimensional_smoke_test
    resultFolder = fullfile(resultFolder, 'smoke');
    figureFolder = fullfile(resultFolder, 'figs');
    fileSuffix = '_smoke';
else
    figureFolder = fullfile(root, 'figs');
    fileSuffix = '';
end
if ~isfolder(resultFolder), mkdir(resultFolder); end
if ~isfolder(figureFolder), mkdir(figureFolder); end

common = struct('grid', grid, 'plan', plan, 'beta', beta, ...
    'delta', delta, 'mass', mass, 'fisher_regularization', fisher, ...
    'potential_regularization', potential, ...
    'regularization', struct('name', 'shift_smooth', ...
        'epsilon', epsilon, 'transition_width', epsilon), ...
    'entropy', struct('enabled', false, 'eta', 0));
parameters = struct('beta', beta, 'delta', delta, 'mass', mass, ...
    'epsilon', epsilon, 'sigma', sigma, 'L', L, 'Nx', Nx, 'Ny', Ny, ...
    'gamma_x', gamma_x, 'gamma_y', gamma_y, 'V0', V0, ...
    'k', lattice_k, 'R0', R0, 'R1', R1, ...
    'R0_fraction', R0_fraction, 'R1_fraction', R1_fraction, ...
    'optical_period_count', n_period, ...
    'smoke_test', two_dimensional_smoke_test);

cases(1).name = 'anisotropic_harmonic';
cases(1).label = 'Anisotropic harmonic trap';
cases(1).V = V_ah;
cases(1).mat_file = fullfile(resultFolder, ...
    ['anisotropic_harmonic_2d', fileSuffix, '.mat']);
cases(1).fig_file = fullfile(figureFolder, ...
    ['ground_state_2d_anisotropic_harmonic', fileSuffix, '.fig']);
cases(1).eps_file = fullfile(figureFolder, ...
    ['ground_state_2d_anisotropic_harmonic', fileSuffix, '.eps']);
cases(2).name = 'optical_lattice';
cases(2).label = 'Harmonic trap with optical lattice';
cases(2).V = V_optical;
cases(2).mat_file = fullfile(resultFolder, ...
    ['optical_lattice_2d', fileSuffix, '.mat']);
cases(2).fig_file = fullfile(figureFolder, ...
    ['ground_state_2d_optical_lattice', fileSuffix, '.fig']);
cases(2).eps_file = fullfile(figureFolder, ...
    ['ground_state_2d_optical_lattice', fileSuffix, '.eps']);

solutions = cell(2, 1);
archives = cell(2, 1);
for caseIndex = 1:2
    fprintf('\n2D case %d/2: %s, %d x %d\n', ...
        caseIndex, cases(caseIndex).label, Nx, Ny);
    compatible = false;
    if reuse_two_dimensional_results && isfile(cases(caseIndex).mat_file)
        saved = load(cases(caseIndex).mat_file, 'archive');
        if isfield(saved, 'archive')
            compatible = compatibleArchive(saved.archive, parameters, ...
                cases(caseIndex).name);
        end
    end
    if compatible
        fprintf('Reusing compatible completed result: %s\n', ...
            cases(caseIndex).mat_file);
        archive = saved.archive;
        rho = archive.rho(:);
    else
        problem = common;
        problem.V = cases(caseIndex).V(:);
        problem.trapping_potential.label = cases(caseIndex).label;
        solveTimer = tic;
        result = src.SolveGroundState1D(problem, rho0, solver);
        measuredTime = toc(solveTimer);
        rho = result.rho;
        diagnostics = twoDimensionalDiagnostics(reshape(rho, grid.shape), ...
            result, grid, plan, R0);
        diagnostics.measured_solver_time = measuredTime;
        archive = makeArchive(cases(caseIndex), parameters, grid, ...
            rho, result, diagnostics, fisher, potential, continuationX, ...
            continuationY);
        save(cases(caseIndex).mat_file, 'archive', '-v7.3');
    end
    solutions{caseIndex} = reshape(rho, grid.shape);
    archives{caseIndex} = archive;
    if ~two_dimensional_smoke_test
        plotDensity(grid, solutions{caseIndex}, cases(caseIndex).label, ...
            cases(caseIndex).fig_file, cases(caseIndex).eps_file);
    end
    printCaseSummary(cases(caseIndex).label, archive.diagnostics);
end

latticeStateDifference = sqrt(grid.h * sum( ...
    (solutions{2}(:) - solutions{1}(:)) .^ 2));
latticeStateMaxDifference = max(abs(solutions{2}(:) - solutions{1}(:)));
fprintf('\nOptical-lattice modulation diagnostics\n');
fprintf('  L2 state difference  : %.6e\n', latticeStateDifference);
fprintf('  max state difference : %.6e\n', latticeStateMaxDifference);

if ~two_dimensional_smoke_test
    for caseIndex = 1:2
        archive = archives{caseIndex};
        archive.lattice_state_difference = latticeStateDifference;
        archive.lattice_state_max_difference = latticeStateMaxDifference;
        save(cases(caseIndex).mat_file, 'archive', '-v7.3');
    end
end

fprintf('\nTeX-ready parameter summary\n');
fprintf(['  L=%g, Nx=Ny=%d, beta=%g, delta=%g, epsilon=%.1e, ' ...
    'sigma=%.1e\n'], L, Nx, beta, delta, epsilon, sigma);
fprintf('  gamma_x=%g, gamma_y=%g, V0=%g, k=pi/2\n', ...
    gamma_x, gamma_y, V0);
fprintf('  harmonic core error: %.3e\n', corePotentialError);
if ~two_dimensional_smoke_test
    assertFiles(cases);
    manuscriptFigureFolder = fullfile(fileparts(root), 'manuscript', 'figs');
    if isfolder(manuscriptFigureFolder)
        for caseIndex = 1:2
            copyfile(cases(caseIndex).eps_file, manuscriptFigureFolder);
        end
    end
end

if do_resolution_check && ~two_dimensional_smoke_test
    warning(['Optional 1024^2 certification is deliberately not automatic ' ...
        'in this display experiment. Set up a separate confirmed run.']);
end

function solver = productionSolver()
config = experiments.DefaultConfig();
solver = config.solver;
solver.name = 'fista_cd';
solver.splitting = 'potential_prox';
solver.projection_name = 'semismooth';
solver.display = false;
solver.pg_tol = 1e-5;
solver.final_pg_tol = 1e-12;
solver.certification_tol = 1e-9;
solver.max_iter = 200000;
solver.residual_check_interval = 1;
solver.switch.enabled = true;
solver.switch.energy_window = 50;
solver.switch.energy_tol = 1e-12;
solver.switch.consecutive_windows = 2;
solver.switch.min_iter = 200;
solver.switch.pg_entry_tol = 1e-5;
solver.switch.max_main_iter = 20000;
solver.switch.forced_pg_tol = 1e-4;
solver.switch.stop_at_energy_handoff = true;
solver.polish_mode = 'if_needed';
solver.polish.pg_tol = 1e-12;
solver.polish.max_iter = 20;
solver.polish.linear_solver = 'interior_pcg_schur';
solver.polish.preconditioner = 'fd_variable';
solver.polish.allow_pdas_fallback = true;
% Use the same certified projection tolerance in the final active-set
% diagnostic; the tighter helper fallback is unnecessarily fragile when
% summing Nx*Ny entries.
solver.polish.projection_tol = solver.projection_tol;
end

function diagnostics = twoDimensionalDiagnostics(rho, result, grid, plan, R0)
continuationMask = abs(grid.X) >= R0 | abs(grid.Y) >= R0;
diagnostics = result.diagnostics;
diagnostics.continuation_mass = grid.h * sum(rho(continuationMask));
diagnostics.continuation_rho_max = max(rho(continuationMask));
rhoHat = fft2(rho);
kx = abs(plan.integer_modes_x);
ky = abs(plan.integer_modes_y);
tailMask = reshape(ky >= 0.8 * max(ky), [], 1) ...
    | reshape(kx >= 0.8 * max(kx), 1, []);
diagnostics.fft_tail_l2 = norm(rhoHat(tailMask)) / norm(rhoHat(:));
reflectX = rho(:, [1, grid.Nx:-1:2]);
reflectY = rho([1, grid.Ny:-1:2], :);
rhoNorm = sqrt(grid.h * sum(rho(:) .^ 2));
diagnostics.x_reflection_defect = sqrt(grid.h ...
    * sum((rho(:) - reflectX(:)) .^ 2)) / rhoNorm;
diagnostics.y_reflection_defect = sqrt(grid.h ...
    * sum((rho(:) - reflectY(:)) .^ 2)) / rhoNorm;
diagnostics.max_density = max(rho(:));
diagnostics.min_density = min(rho(:));
diagnostics.exact_zero_count = nnz(rho == 0);
if diagnostics.fft_tail_l2 > 1e-6
    warning('N=%d may not sufficiently resolve the 2D state (FFT tail %.3e).', ...
        grid.Nx, diagnostics.fft_tail_l2);
end
if diagnostics.continuation_mass > 1e-8 ...
        || diagnostics.continuation_rho_max > 1e-8
    warning(['Density is not negligible in the potential continuation ' ...
        'region: mass %.3e, max %.3e.'], ...
        diagnostics.continuation_mass, diagnostics.continuation_rho_max);
end
end

function archive = makeArchive(caseData, parameters, grid, rho, result, ...
        diagnostics, fisher, potential, continuationX, continuationY)
archive.case_name = caseData.name;
archive.case_label = caseData.label;
archive.parameters = parameters;
archive.x = grid.x;
archive.y = grid.y;
archive.rho = reshape(rho, grid.shape);
archive.V = caseData.V;
archive.energy = result.target_energy;
archive.PG = diagnostics.final_pg_residual;
archive.KKT = diagnostics.final_kkt_residual;
archive.mass_error = diagnostics.mass_error;
archive.min_rho = diagnostics.min_density;
archive.max_rho = diagnostics.max_density;
archive.continuation_mass = diagnostics.continuation_mass;
archive.continuation_rho_max = diagnostics.continuation_rho_max;
archive.fft_tail_l2 = diagnostics.fft_tail_l2;
archive.diagnostics = diagnostics;
archive.solver_history = result.history;
archive.solver = result.solver;
archive.fisher_regularization = fisher;
archive.potential_regularization = potential;
archive.continuation_x = continuationX;
archive.continuation_y = continuationY;
archive.gamma_x = parameters.gamma_x;
archive.gamma_y = parameters.gamma_y;
if strcmp(caseData.name, 'optical_lattice')
    archive.V0 = parameters.V0;
    archive.k = parameters.k;
end
end

function compatible = compatibleArchive(archive, parameters, caseName)
compatible = isfield(archive, 'parameters') ...
    && isfield(archive, 'case_name') ...
    && strcmp(archive.case_name, caseName);
if ~compatible, return; end
fields = {'beta', 'delta', 'mass', 'epsilon', 'sigma', 'L', 'Nx', ...
    'Ny', 'gamma_x', 'gamma_y', 'V0', 'k', 'R0', 'R1'};
for j = 1:numel(fields)
    name = fields{j};
    compatible = compatible && isfield(archive.parameters, name) ...
        && archive.parameters.(name) == parameters.(name);
end
end

function plotDensity(grid, rho, titleText, figFile, epsFile)
f = figure('Color', 'w', 'Position', [100, 100, 650, 560]);
ax = axes(f);
contourf(ax, grid.X, grid.Y, rho, 32, 'LineStyle', 'none');
axis(ax, 'equal'); axis(ax, 'tight'); box(ax, 'on');
xlim(ax, [-5, 5]);
ylim(ax, [-5, 5]);
colormap(ax, jet(256));
cb = colorbar(ax);
set(cb, 'FontSize', 20);
set(ax, 'FontSize', 20, 'LineWidth', 1.2);
xlabel(ax, '$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel(ax, '$y$', 'Interpreter', 'latex', 'FontSize', 20);
title(ax, titleText, 'Interpreter', 'latex', ...
    'FontSize', 20, 'FontWeight', 'normal');
savefig(f, figFile);
screenDpi = get(groot, 'ScreenPixelsPerInch');
figurePosition = get(f, 'Position');
paperSize = figurePosition(3:4) / screenDpi;
set(f, 'Renderer', 'painters', 'PaperUnits', 'inches', ...
    'PaperSize', paperSize, 'PaperPosition', [0, 0, paperSize], ...
    'PaperPositionMode', 'manual');
print(f, epsFile, '-depsc2', '-painters');
close(f);
end

function printCaseSummary(label, diagnostics)
fprintf('%s\n', label);
fprintf('  energy                 : %.15e\n', diagnostics.final_energy);
fprintf('  PG / KKT               : %.3e / %.3e\n', ...
    diagnostics.final_pg_residual, diagnostics.final_kkt_residual);
fprintf('  mass error             : %.3e\n', diagnostics.mass_error);
fprintf('  min / max rho          : %.3e / %.3e\n', ...
    diagnostics.min_density, diagnostics.max_density);
fprintf('  exact zeros            : %d\n', diagnostics.exact_zero_count);
fprintf('  FISTA / Newton         : %d / %d\n', ...
    diagnostics.main_iterations, diagnostics.polish_iterations);
fprintf('  total PCG              : %d\n', ...
    diagnostics.total_pcg_z_iterations ...
    + diagnostics.total_pcg_w_iterations);
fprintf('  solver time            : %.3f s\n', ...
    diagnostics.measured_solver_time);
fprintf('  continuation mass/max  : %.3e / %.3e\n', ...
    diagnostics.continuation_mass, diagnostics.continuation_rho_max);
fprintf('  FFT tail               : %.3e\n', diagnostics.fft_tail_l2);
fprintf('  reflection defects x/y : %.3e / %.3e\n', ...
    diagnostics.x_reflection_defect, diagnostics.y_reflection_defect);
end

function assertFiles(cases)
for j = 1:numel(cases)
    names = {cases(j).mat_file, cases(j).fig_file, cases(j).eps_file};
    for k = 1:numel(names)
        info = dir(names{k});
        assert(~isempty(info) && info.bytes > 0, ...
            'Required output is missing or empty: %s', names{k});
    end
end
end
