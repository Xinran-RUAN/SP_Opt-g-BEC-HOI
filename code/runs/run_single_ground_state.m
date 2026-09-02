%RUN_SINGLE_GROUND_STATE Solve one configurable discrete ground state.
clearvars; clc;

% ======================== regularization ===============================
epsilon = 1e-3;

% Fisher: s_epsilon(rho)
s_epsilon = @(rho) rho + epsilon;
ds_epsilon = @(rho) ones(size(rho));
d2s_epsilon = @(rho) zeros(size(rho));
fisher_label = 's_epsilon(rho) = rho + epsilon';

% Potential: p_sigma(rho). Set sigma directly.
sigma = 1e-12;
p_sigma = @(rho) rho .^ 2 ./ (hypot(rho, sigma) + sigma);
dp_sigma = @(rho) rho ./ hypot(rho, sigma);
d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
p_sigma_label = 'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
% =======================================================================

beta = 10;
delta = 10;
L = 32;
N = 128;

% ===================== trapping potential V(x) =========================
trapping_potential_choice = 'harmonic';
% Set to 'harmonic_cinf_periodic' for an unchanged harmonic core and a
% boundary-only C-infinity periodic continuation.
R0_fraction = 0.75;
R1_fraction = 0.90;
negligible_density_tol = 1e-10; % diagnostic only
negligible_mass_tol = 1e-10;    % diagnostic only
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

entropy_enabled = false;
eta = 0;

solver_name = 'fista_cd';
projection_name = 'semismooth';

pg_tol = 1e-8;
max_iter = 200000;
polish_mode = 'if_needed';
polish_pg_tol = 1e-12;

show_plot = true;
save_result = true;
overwrite_existing = true;
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_HOI();

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
config.parameters.L = L;
config.parameters.N = N;
config.parameters.epsilon = epsilon;
config.fisher_regularization.epsilon = epsilon;
config.fisher_regularization.s_epsilon = s_epsilon;
config.fisher_regularization.ds_epsilon = ds_epsilon;
config.fisher_regularization.d2s_epsilon = d2s_epsilon;
config.fisher_regularization.label = fisher_label;
config.potential_regularization.sigma = sigma;
config.potential_regularization.name = 'inline_p_sigma'; % metadata only
config.potential_regularization.p_sigma = p_sigma;
config.potential_regularization.dp_sigma = dp_sigma;
config.potential_regularization.d2p_sigma = d2p_sigma;
config.potential_regularization.label = p_sigma_label;
config.potential_regularization.prox_type = 'generic_convex';
config.entropy.enabled = entropy_enabled;
config.entropy.eta = eta;
if strcmpi(solver_name, 'auto')
    if entropy_enabled && eta > 0
        solver_name = 'fista_cd';
    else
        solver_name = 'spg';
    end
end
config.solver.name = solver_name;
config.solver.splitting = 'potential_prox';
config.solver.projection_name = projection_name;
config.solver.pg_tol = pg_tol;
config.solver.max_iter = max_iter;
config.solver.polish_mode = polish_mode;
config.solver.polish.pg_tol = polish_pg_tol;
config.solver.switch.enabled = true;
config.output.show_plot = show_plot;
config.output.save_result = save_result;
config.output.overwrite_existing = overwrite_existing;

result = experiments.SolveGroundState(config);

active = result.diagnostics.active_set;
fourier = result.diagnostics.fourier_tail;
[~, nearIndex] = min(abs(active.zero_tol_list - 1e-12));
fprintf('\nRegularized density energy E_{epsilon,sigma}\n');
fprintf('V(x)           : %s\n', V_label);
fprintf('epsilon        : %.3e\n', epsilon);
fprintf('s_epsilon      : %s\n', fisher_label);
fprintf('sigma          : %.3e\n', sigma);
fprintf('p_sigma        : %s\n', p_sigma_label);
fprintf('Entropy        : %s\n', onOff(entropy_enabled && eta > 0));
fprintf('eta            : %.3e\n', eta);
fprintf('Main solver    : %s\n', result.diagnostics.main_solver);
fprintf('Main iter      : %d\n', result.diagnostics.main_iterations);
fprintf('Target E_eps,sigma : %.15e\n', result.target_energy);
fprintf('Baseline E_eps,0   : %.15e\n', result.baseline_energy);
fprintf('Stationarity   : %.3e\n', result.diagnostics.pg_residual);
fprintf('KKT residual   : %.3e\n', result.diagnostics.kkt_residual);
fprintf('Min density    : %.3e\n', min(result.rho));
fprintf('Exact zeros    : %d\n', active.exact_zero_count);
fprintf('Near-zero frac : %.3e\n', active.zero_fraction_by_tol(nearIndex));
fprintf('FFT tail ratio : %.3e\n', fourier.tail_ratio_quarter);
fprintf('Polish         : %s\n', result.diagnostics.polish_status);
fprintf('Polish iter    : %d\n', result.diagnostics.polish_iterations);
fprintf('Mass error     : %.3e\n', result.diagnostics.mass_error);
fprintf('Tail mass      : %.3e\n', result.diagnostics.tail_mass);
modification = result.diagnostics.trapping_potential_modification;
if modification.enabled
    fprintf('Harmonic core error : %.3e\n', ...
        modification.core_potential_error);
    fprintf('Modified rho max    : %.3e\n', modification.rho_max);
    fprintf('Modified mass       : %.3e\n', modification.mass);
    if ~modification.density_negligible ...
            || ~modification.mass_negligible
        warning(['The boundary continuation acts where density is not ' ...
            'negligible at the configured diagnostic tolerances.']);
    end
end

if show_plot
    figure('Name', '1D HOI ground-state density');
    plot(result.grid.x, result.rho, 'LineWidth', 1.4);
    xlabel('x'); ylabel('\rho'); grid on;
    title(sprintf('E_{eps,sigma}, N=%d, epsilon=%.1e, sigma=%.1e', ...
        N, epsilon, sigma), ...
        'Interpreter', 'none');
end

if save_result
    baseName = sprintf('single_inline_sigma_%s_N%d_eps_%s.mat', ...
        numberTag(sigma), N, strrep(sprintf('%.0e', epsilon), '-', 'm'));
    savedFile = experiments.SaveResult( ...
        result, 'single', baseName, overwrite_existing);
    if ~isempty(savedFile)
        fprintf('saved: %s\n', savedFile);
    end
end

function tag = numberTag(value)
tag = strrep(strrep(sprintf('%.0e', value), '-', 'm'), '+', 'p');
end

function label = onOff(active)
if active
    label = 'on';
else
    label = 'off';
end
end
