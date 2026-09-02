function result = SolveGroundState1D(problem, rho0, solver)
%SOLVEGROUNDSTATE1D Convex minimizer followed by optional KKT polish.

requiredProblemFields = {'grid', 'plan', 'V', 'beta', 'delta', ...
    'mass'};
if ~isstruct(problem) || ~all(isfield(problem, requiredProblemFields))
    error('src:SolveGroundState1D:InvalidProblem', ...
        'problem is missing required fields.');
end
if numel(problem.V) ~= problem.grid.N || problem.plan.N ~= problem.grid.N
    error('src:SolveGroundState1D:SizeMismatch', ...
        'Grid, plan, and potential sizes must agree.');
end
[problem, ~] = applyFisherRegularization(problem);

problem = applyEntropyDefaults(problem);
[problem, ~] = ...
    applyPotentialRegularization(problem);
solver = applySolverDefaults(solver);
entropyActive = problem.entropy.enabled && problem.entropy.eta > 0;
if entropyActive && strcmpi(solver.name, 'spg')
    error('src:SolveGroundState1D:SPGNotAvailableForEntropy', ...
        'Classic SPG is retained only for eta=0; use fista_cd or ista.');
end
if strcmpi(solver.splitting, 'potential_prox') ...
        && ~strcmpi(solver.name, 'fista_cd') ...
        && ~strcmpi(solver.name, 'fista-cd')
    error('src:SolveGroundState1D:PotentialProxRequiresFISTA', ...
        'potential_prox splitting is implemented for FISTA-CD only.');
end
rho0 = projectInitial(rho0, problem, solver);
[problem.fisher_regularization, fisherRegularizationInfo] = ...
    src.regularization.ValidateFisher( ...
    problem.fisher_regularization, [0; rho0(:)]);
[problem.potential_regularization, potentialRegularizationInfo] = ...
    src.potential.Validate(problem.potential_regularization, ...
    problem.fisher_regularization.epsilon, [0; rho0(:)]);
printMainHeader(solver);

switch lower(char(solver.name))
    case 'ista'
        [rhoMain, mainDiagnostics, mainHistory] = ...
            src.solvers.ISTA(problem, rho0, solver);
    case 'spg'
        [rhoMain, mainDiagnostics, mainHistory] = ...
            src.solvers.SPG(problem, rho0, solver);
    case {'fista_cd', 'fista-cd'}
        solver.name = 'fista_cd';
        [rhoMain, mainDiagnostics, mainHistory] = ...
            src.solvers.FISTACD(problem, rho0, solver);
    otherwise
        error('src:SolveGroundState1D:UnknownSolver', ...
            'Unknown solver "%s".', solver.name);
end

mainState = src.solvers.EvaluateState(rhoMain, problem, solver);
mainFullPg = fullResidual(rhoMain, mainState, problem, solver);
polishResult = emptyPolishResult( ...
    rhoMain, mainState, mainFullPg, solver.polish.pg_tol);
switch lower(char(solver.polish_mode))
    case 'none'
        shouldPolish = false;
    case 'always'
        shouldPolish = true;
    case 'if_needed'
        if isfield(mainDiagnostics, 'request_polish') ...
                && solver.switch.enabled
            shouldPolish = mainDiagnostics.request_polish;
        else
            shouldPolish = mainFullPg > solver.polish.pg_tol;
        end
    otherwise
        error('src:SolveGroundState1D:UnknownPolishMode', ...
            'polish_mode must be none, always, or if_needed.');
end

if shouldPolish && ~mainDiagnostics.failed
    polishOptions = solver.polish;
    polishOptions.projection_name = solver.projection_name;
    polishOptions.projection_tol = solver.projection_tol;
    polishOptions.residual_step = solver.residual_step;
    polishOptions.display = solver.display;
    if solver.switch.enabled
        polishOptions.entry_pg_tol = max( ...
            solver.switch.pg_entry_tol, solver.switch.forced_pg_tol);
    end
    if solver.display && solver.switch.enabled
        fprintf('\n[FISTA -> KKT POLISH]\n');
        fprintf('iter       = %d\n', mainDiagnostics.iterations);
        fprintf('energy     = %.15e\n', mainState.energy);
        fprintf('PG         = %.3e\n', mainFullPg);
        fprintf('window dE  = %.3e\n', getFieldOr( ...
            mainDiagnostics, 'main_energy_window_span', NaN));
        fprintf('reason     = %s\n\n', getFieldOr( ...
            mainDiagnostics, 'main_stop_reason', 'polish_requested'));
    end
    if entropyActive
        if solver.display
            fprintf('Entropy Newton   pg_res      kkt_res     GMRES  step\n');
        end
        polishResult = src.solvers.PolishEntropyNewton( ...
            rhoMain, problem, polishOptions);
    else
        polishResult = src.solvers.PolishKKT( ...
            rhoMain, problem, polishOptions);
    end
end

rhoFinal = polishResult.rho;
finalState = src.solvers.EvaluateState(rhoFinal, problem, solver);
diagnostics = mainDiagnostics;
diagnostics.main_solver = upper(strrep(solver.name, '_', '-'));
diagnostics.main_iterations = mainDiagnostics.iterations;
diagnostics.main_elapsed_time = mainDiagnostics.elapsed_time;
diagnostics.main_converged = mainDiagnostics.converged;
diagnostics.main_stop_reason = getFieldOr( ...
    mainDiagnostics, 'main_stop_reason', 'solver_stop');
diagnostics.main_energy = mainState.energy;
diagnostics.main_pg_residual = mainFullPg;
diagnostics.main_energy_window_span = getFieldOr( ...
    mainDiagnostics, 'main_energy_window_span', NaN);
diagnostics.energy_before_polish = mainState.energy;
diagnostics.physical_energy_before_polish = mainState.physical_energy;
diagnostics.entropy_value_before_polish = mainState.entropy_value;
diagnostics.augmented_energy_before_polish = mainState.augmented_energy;
diagnostics.pg_before_polish = mainFullPg;
diagnostics.kkt_before_polish = mainState.kkt_residual;
diagnostics.polish_mode = solver.polish_mode;
diagnostics.polish_entered = polishResult.polish_attempted;
diagnostics.polish_attempted = polishResult.polish_attempted;
diagnostics.polish_converged = polishResult.polish_converged;
diagnostics.polish_failed = polishResult.failed;
diagnostics.polish_failure_message = polishResult.failure_message;
diagnostics.polish_status = polishResult.status;
diagnostics.polish_iterations = polishResult.iterations;
diagnostics.polish_elapsed_time = polishResult.elapsed_time;
diagnostics.polish_solver = getFieldOr( ...
    polishResult, 'polish_solver', 'none');
diagnostics.polish_preconditioner = getFieldOr( ...
    polishResult, 'preconditioner', 'none');
diagnostics.polish_interior_residual = getFieldOr( ...
    polishResult, 'interior_residual', NaN);
diagnostics.total_pcg_z_iterations = getFieldOr( ...
    polishResult, 'total_pcg_z_iterations', 0);
diagnostics.total_pcg_w_iterations = getFieldOr( ...
    polishResult, 'total_pcg_w_iterations', 0);
diagnostics.max_pcg_iterations = getFieldOr( ...
    polishResult, 'max_pcg_iterations', 0);
diagnostics.total_gmres_iterations = getFieldOr( ...
    polishResult, 'total_gmres_iterations', 0);
diagnostics.max_gmres_iterations = getFieldOr( ...
    polishResult, 'max_gmres_iterations', 0);
diagnostics.polish_switched_to_pdas = getFieldOr( ...
    polishResult, 'switched_to_pdas', false);
diagnostics.polish_fallback_reason = getFieldOr( ...
    polishResult, 'fallback_reason', '');
diagnostics.projected_active_count_handoff = getFieldOr( ...
    polishResult, 'handoff_projected_active_count', NaN);
diagnostics.energy_after_polish = finalState.energy;
diagnostics.physical_energy_after_polish = finalState.physical_energy;
diagnostics.entropy_value_after_polish = finalState.entropy_value;
diagnostics.augmented_energy_after_polish = finalState.augmented_energy;
diagnostics.pg_after_polish = finalState.pg_residual;
diagnostics.kkt_after_polish = finalState.kkt_residual;
diagnostics.polish_energy_change = finalState.energy - mainState.energy;
diagnostics.polish_state_change = sqrt(problem.grid.h ...
    * sum((rhoFinal - rhoMain) .^ 2));
diagnostics.energy = finalState.energy;
diagnostics.target_energy = finalState.energy;
diagnostics.baseline_energy = src.discretization.ps.BaselineEnergy( ...
    rhoFinal, problem);
diagnostics.potential_type = problem.potential_regularization.name;
diagnostics.potential_power = problem.potential_regularization.power;
diagnostics.potential_sigma = problem.potential_regularization.sigma;
diagnostics.potential_label = problem.potential_regularization.label;
if isfield(problem, 'trapping_potential') ...
        && isfield(problem.trapping_potential, 'label')
    diagnostics.trapping_potential_label = ...
        problem.trapping_potential.label;
elseif isfield(problem, 'potential_label')
    diagnostics.trapping_potential_label = problem.potential_label;
else
    diagnostics.trapping_potential_label = 'nodal V';
end
diagnostics.fisher_label = problem.fisher_regularization.label;
diagnostics.fisher_epsilon = problem.fisher_regularization.epsilon;
diagnostics.physical_energy = finalState.physical_energy;
diagnostics.entropy_value = finalState.entropy_value;
diagnostics.augmented_energy = finalState.augmented_energy;
diagnostics.entropy_enabled = entropyActive;
diagnostics.entropy_eta = problem.entropy.eta;
diagnostics.pg_residual = finalState.pg_residual;
if entropyActive
    diagnostics.composite_pg_residual = finalState.pg_residual;
    diagnostics.full_pg_residual = NaN;
else
    diagnostics.composite_pg_residual = ...
        src.solvers.CompositeGradientMapping( ...
        rhoFinal, problem, solver);
    diagnostics.full_pg_residual = src.solvers.FullGradientMapping( ...
        rhoFinal, finalState.gradient, problem, solver);
end
diagnostics.final_energy = finalState.energy;
diagnostics.final_pg_residual = fullResidual( ...
    rhoFinal, finalState, problem, solver);
diagnostics.final_kkt_residual = finalState.kkt_residual;
if ~entropyActive
    diagnostics.pg_residual = diagnostics.final_pg_residual;
end
diagnostics.kkt_residual = finalState.kkt_residual;
diagnostics.exact_zero_count = nnz(rhoFinal == 0);
diagnostics.strict_positive_count = nnz(rhoFinal > 0);
if ~entropyActive
    projectedDiagnostic = src.solvers.InteriorKKTResidual( ...
        rhoFinal, problem, solver.polish, finalState.gradient);
    diagnostics.projected_active_count_final = ...
        projectedDiagnostic.projected_active_count;
    diagnostics.projected_free_count_final = ...
        projectedDiagnostic.projected_free_count;
else
    diagnostics.projected_active_count_final = NaN;
    diagnostics.projected_free_count_final = NaN;
end
diagnostics.mass_error = finalState.mass_error;
diagnostics.min_density = finalState.min_density;
diagnostics.total_elapsed_time = mainDiagnostics.elapsed_time ...
    + polishResult.elapsed_time;
diagnostics.total_iterations = mainDiagnostics.iterations ...
    + polishResult.iterations;
if solver.switch.enabled
    finalTolerance = solver.final_pg_tol;
elseif strcmpi(solver.polish_mode, 'none')
    finalTolerance = solver.pg_tol;
else
    finalTolerance = solver.polish.pg_tol;
end
diagnostics.converged = diagnostics.final_pg_residual <= finalTolerance;
if entropyActive
    diagnostics.certified = diagnostics.converged;
else
    diagnostics.certified = diagnostics.full_pg_residual ...
        <= solver.certification_tol;
end
if ~solver.switch.enabled && strcmpi(solver.splitting, 'potential_prox') ...
        && diagnostics.composite_pg_residual <= solver.pg_tol ...
        && diagnostics.full_pg_residual > solver.certification_tol
    warning('src:SolveGroundState1D:FullResidualNotCertified', ...
        ['Composite PG converged, but full PG %.3e exceeds the ' ...
        'certification tolerance %.3e.'], ...
        diagnostics.full_pg_residual, solver.certification_tol);
end

result.rho = rhoFinal;
result.energy = finalState.energy;
result.target_energy = finalState.energy;
result.baseline_energy = diagnostics.baseline_energy;
result.physical_energy = finalState.physical_energy;
result.entropy_value = finalState.entropy_value;
result.augmented_energy = finalState.augmented_energy;
result.diagnostics = diagnostics;
result.history.main = mainHistory;
result.history.polish = polishResult.history;
result.solver = solver;
result.regularization = fisherRegularizationInfo;
result.fisher_regularization = fisherRegularizationInfo;
result.problem = problem;
result.main_rho = rhoMain;
result.entropy = problem.entropy;
result.potential_regularization = potentialRegularizationInfo;

if solver.display
    printSummary(diagnostics);
end
end

function solver = applySolverDefaults(solver)
if nargin < 1 || isempty(solver)
    solver = struct();
end

defaults.name = 'spg';
defaults.splitting = 'legacy_full_gradient';
defaults.projection_name = 'simplex';
defaults.projection_tol = 1e-13;
defaults.residual_step = 1;
defaults.residual_check_interval = 1;
defaults.active_tol = 1e-12;
defaults.L0 = 1;
defaults.backtrack_factor = 2;
defaults.ista_L_decrease = 0.5;
defaults.max_backtracks = 100;
defaults.pg_tol = 1e-8;
defaults.certification_tol = 1e-9;
defaults.final_pg_tol = 1e-12;
defaults.fallback_fista_iter = 0;
defaults.energy_tol = 1e-14;
defaults.max_iter = 200000;
defaults.a = 4;
defaults.feasibility_tol = 1e-13;
defaults.display = true;
defaults.display_interval = 100;
defaults.display_every = 200;
defaults.spg_bb_type = 'bb1';
defaults.spg_alpha_min = 1e-12;
defaults.spg_alpha_max = 1e12;
defaults.spg_alpha_reset = 1;
defaults.spg_curvature_tol = 1e-18;
defaults.spg_nonmonotone_M = 10;
defaults.spg_c1 = 1e-4;
defaults.spg_backtrack = 0.5;
defaults.spg_min_lambda = 1e-14;
defaults.spg_descent_tol = 1e-13;
defaults.polish_mode = 'if_needed';

names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(solver, names{j}) || isempty(solver.(names{j}))
        solver.(names{j}) = defaults.(names{j});
    end
end
if ~isfield(solver, 'switch') || isempty(solver.switch)
    solver.switch = struct();
end
switchDefaults.enabled = false;
switchDefaults.energy_window = 50;
switchDefaults.energy_tol = 1e-12;
switchDefaults.consecutive_windows = 2;
switchDefaults.min_iter = 200;
switchDefaults.pg_entry_tol = 1e-5;
switchDefaults.max_main_iter = 20000;
switchDefaults.forced_pg_tol = 1e-4;
switchDefaults.stop_at_energy_handoff = true;
names = fieldnames(switchDefaults);
for j = 1:numel(names)
    if ~isfield(solver.switch, names{j}) ...
            || isempty(solver.switch.(names{j}))
        solver.switch.(names{j}) = switchDefaults.(names{j});
    end
end
if ~isfield(solver, 'history') || isempty(solver.history)
    solver.history = struct();
end
if ~isfield(solver.history, 'capture_handoff_state') ...
        || isempty(solver.history.capture_handoff_state)
    solver.history.capture_handoff_state = false;
end
if ~isfield(solver, 'potential_prox') || isempty(solver.potential_prox)
    solver.potential_prox = struct();
end
potentialProxDefaults.mass_tol = 1e-13;
potentialProxDefaults.lambda_max_iter = 50;
potentialProxDefaults.inner_tol = 1e-13;
potentialProxDefaults.inner_max_iter = 30;
potentialProxDefaults.bracket_max_iter = 60;
names = fieldnames(potentialProxDefaults);
for j = 1:numel(names)
    if ~isfield(solver.potential_prox, names{j}) ...
            || isempty(solver.potential_prox.(names{j}))
        solver.potential_prox.(names{j}) = potentialProxDefaults.(names{j});
    end
end
if ~isfield(solver, 'polish') || isempty(solver.polish)
    solver.polish = struct();
end
polishDefaults.entry_pg_tol = 1e-5;
polishDefaults.pg_tol = 1e-12;
polishDefaults.max_iter = 20;
polishDefaults.linear_solver = 'interior_pcg_schur';
polishDefaults.preconditioner = 'fd_variable';
polishDefaults.allow_pdas_fallback = true;
polishDefaults.pdas_max_iter = 50;
polishDefaults.tau_active = 1;
polishDefaults.active_step = 1;
polishDefaults.active_tol = 1e-12;
polishDefaults.gmres_tol = 1e-10;
polishDefaults.gmres_maxit = 200;
polishDefaults.gmres_failure_tol = 0.5;
polishDefaults.pcg_tol_max = 1e-2;
polishDefaults.pcg_tol_min = 1e-10;
polishDefaults.pcg_forcing_factor = 0.1;
polishDefaults.pcg_maxit = [];
polishDefaults.fraction_to_boundary = 0.995;
polishDefaults.residual_armijo = 1e-4;
polishDefaults.stagnation_window = 5;
polishDefaults.stagnation_rel_improvement = 1e-2;
polishDefaults.acceptable_floor = 1e-10;
polishDefaults.residual_c1 = 1e-4;
polishDefaults.backtrack = 0.5;
polishDefaults.min_step = 1e-12;
polishDefaults.max_backtracks = 30;
polishDefaults.mass_tol = 1e-12;
polishDefaults.roundoff_tol = 1e-13;
names = fieldnames(polishDefaults);
for j = 1:numel(names)
    if ~isfield(solver.polish, names{j}) || isempty(solver.polish.(names{j}))
        solver.polish.(names{j}) = polishDefaults.(names{j});
    end
end

if ~ismember(lower(char(solver.projection_name)), {'simplex', 'semismooth'})
    error('src:SolveGroundState1D:UnknownProjection', ...
        'projection_name must be simplex or semismooth.');
end
if ~ismember(lower(char(solver.splitting)), ...
        {'potential_prox', 'legacy_full_gradient'})
    error('src:SolveGroundState1D:UnknownSplitting', ...
        'splitting must be potential_prox or legacy_full_gradient.');
end
if ~ismember(lower(char(solver.spg_bb_type)), {'bb1', 'bb2', 'alternate'})
    error('src:SolveGroundState1D:UnknownBBType', ...
        'spg_bb_type must be bb1, bb2, or alternate.');
end
if solver.L0 <= 0 || solver.backtrack_factor <= 1 ...
        || solver.ista_L_decrease <= 0 || solver.ista_L_decrease > 1 ...
        || solver.max_iter < 1 || solver.pg_tol <= 0 ...
        || solver.certification_tol <= 0 || solver.final_pg_tol <= 0 ...
        || solver.fallback_fista_iter < 0 || solver.a <= 2 ...
        || solver.residual_step <= 0 ...
        || solver.residual_check_interval < 1 ...
        || solver.residual_check_interval ...
        ~= round(solver.residual_check_interval) ...
        || solver.spg_alpha_min <= 0 ...
        || solver.spg_alpha_max < solver.spg_alpha_min ...
        || ~isfinite(solver.spg_alpha_reset) || solver.spg_alpha_reset <= 0 ...
        || solver.spg_nonmonotone_M < 1 ...
        || solver.spg_nonmonotone_M ~= round(solver.spg_nonmonotone_M) ...
        || solver.spg_c1 <= 0 || solver.spg_c1 >= 1 ...
        || solver.spg_min_lambda <= 0 ...
        || solver.spg_backtrack <= 0 || solver.spg_backtrack >= 1
    error('src:SolveGroundState1D:InvalidSolverOptions', ...
        'Invalid solver tolerances, residual step, or line-search parameters.');
end
if solver.switch.energy_window < 2 ...
        || solver.switch.energy_window ~= round(solver.switch.energy_window) ...
        || solver.switch.energy_tol <= 0 ...
        || solver.switch.consecutive_windows < 1 ...
        || solver.switch.consecutive_windows ...
        ~= round(solver.switch.consecutive_windows) ...
        || solver.switch.min_iter < 1 ...
        || solver.switch.max_main_iter < solver.switch.min_iter ...
        || solver.switch.pg_entry_tol <= 0 ...
        || solver.switch.forced_pg_tol <= 0 ...
        || ~isLogicalScalar(solver.switch.stop_at_energy_handoff) ...
        || ~isLogicalScalar(solver.history.capture_handoff_state)
    error('src:SolveGroundState1D:InvalidSwitchOptions', ...
        'Invalid FISTA-to-polish switch options.');
end
end

function valid = isLogicalScalar(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) ...
    && isfinite(value) && (value == 0 || value == 1);
end

function problem = applyEntropyDefaults(problem)
if ~isfield(problem, 'entropy') || isempty(problem.entropy)
    problem.entropy = struct();
end

defaults.enabled = false;
defaults.eta = 0;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(problem.entropy, names{j}) ...
            || isempty(problem.entropy.(names{j}))
        problem.entropy.(names{j}) = defaults.(names{j});
    end
end
if ~islogical(problem.entropy.enabled) && ...
        ~(isnumeric(problem.entropy.enabled) && isscalar(problem.entropy.enabled))
    error('src:SolveGroundState1D:InvalidEntropyEnabled', ...
        'entropy.enabled must be a logical scalar.');
end
if ~isscalar(problem.entropy.eta) || ~isfinite(problem.entropy.eta) ...
        || problem.entropy.eta < 0
    error('src:SolveGroundState1D:InvalidEntropyEta', ...
        'entropy.eta must be a finite nonnegative scalar.');
end
if ~isfield(problem.entropy, 'prox') || isempty(problem.entropy.prox)
    problem.entropy.prox = struct();
end
proxDefaults.mass_tol = 1e-14;
proxDefaults.lambda_max_iter = 50;
proxDefaults.inner_newton_max_iter = 30;
proxDefaults.inner_tol = 1e-14;
names = fieldnames(proxDefaults);
for j = 1:numel(names)
    if ~isfield(problem.entropy.prox, names{j}) ...
            || isempty(problem.entropy.prox.(names{j}))
        problem.entropy.prox.(names{j}) = proxDefaults.(names{j});
    end
end
end

function [problem, info] = applyPotentialRegularization(problem)
if ~isfield(problem, 'potential_regularization') ...
        || isempty(problem.potential_regularization)
    problem.potential_regularization.name = 'linear';
end
if ~isfield(problem.potential_regularization, 'convexity_tol') ...
        || isempty(problem.potential_regularization.convexity_tol)
    problem.potential_regularization.convexity_tol = 1e-14;
end
convexityTolerance = problem.potential_regularization.convexity_tol;
if ~isscalar(convexityTolerance) || ~isfinite(convexityTolerance) ...
        || convexityTolerance < 0
    error('src:SolveGroundState1D:InvalidPotentialConvexityTolerance', ...
        'potential_regularization.convexity_tol must be nonnegative.');
end
[problem.potential_regularization, info] = src.potential.Validate( ...
    problem.potential_regularization, ...
    problem.fisher_regularization.epsilon);
if ~strcmpi(problem.potential_regularization.prox_type, 'linear') ...
        && min(problem.V(:)) < -convexityTolerance
    error('src:SolveGroundState1D:NegativePotentialRegularization', ...
        ['Potential regularization requires nonnegative V ' ...
        'to preserve convexity.']);
end
info.convexity_tol = convexityTolerance;
info.min_potential = min(problem.V(:));
end

function [problem, info] = applyFisherRegularization(problem)
if isfield(problem, 'fisher_regularization') ...
        && ~isempty(problem.fisher_regularization)
    fisher = problem.fisher_regularization;
elseif isfield(problem, 'regularization')
    [problem.regularization, ~] = ...
        src.regularization.Validate(problem.regularization);
    fisher = src.regularization.MakeBuiltIn(problem.regularization);
else
    error('src:SolveGroundState1D:MissingFisherRegularization', ...
        ['problem must provide fisher_regularization handles or a ' ...
        'legacy regularization configuration.']);
end
[fisher, info] = src.regularization.ValidateFisher(fisher, ...
    [0; problem.mass / problem.grid.domain_length; ...
    problem.mass / problem.grid.h]);
problem.fisher_regularization = fisher;
if ~isfield(problem, 'regularization')
    problem.regularization.name = 'inline_handles';
    problem.regularization.epsilon = fisher.epsilon;
    problem.regularization.transition_width = fisher.epsilon;
end
end

function rho = projectInitial(rho0, problem, solver)
if numel(rho0) ~= problem.grid.N || any(~isfinite(rho0(:)))
    error('src:SolveGroundState1D:InvalidInitialDensity', ...
        'Initial density must be finite and have grid.N entries.');
end
switch lower(char(solver.projection_name))
    case 'simplex'
        rho = src.constraints.ProjectSimplex( ...
            rho0, problem.mass, problem.grid.h);
    case 'semismooth'
        rho = src.constraints.ProjectPositiveConservative( ...
            rho0, problem.mass, problem.grid.h, solver.projection_tol);
end
end

function polish = emptyPolishResult(rho, state, fullPg, requestedTolerance)
polish.rho = rho;
polish.energy = state.energy;
polish.history = struct();
polish.polish_attempted = false;
polish.polish_converged = fullPg <= requestedTolerance;
polish.failed = false;
polish.failure_message = '';
polish.iterations = 0;
polish.elapsed_time = 0;
polish.status = 'not_attempted';
end

function residual = fullResidual(rho, state, problem, solver)
if problem.entropy.enabled && problem.entropy.eta > 0
    residual = state.pg_residual;
else
    residual = src.solvers.FullGradientMapping( ...
        rho, state.gradient, problem, solver);
end
end

function value = getFieldOr(record, name, fallback)
if isfield(record, name)
    value = record.(name);
else
    value = fallback;
end
end

function printMainHeader(solver)
if ~solver.display
    return;
end
if strcmpi(solver.name, 'spg')
    fprintf('   iter   energy             pg_res       dE          alphaBB     lambda      bt\n');
else
    if solver.switch.enabled
        fprintf('   iter   energy             full_PG      win_dE      L           restart\n');
    else
        fprintf('   iter   energy             pg_res       dE          L           restart\n');
    end
end
end

function printSummary(diagnostics)
fprintf('\nMain solver : %s\n', diagnostics.main_solver);
fprintf('Main iter   : %d\n', diagnostics.main_iterations);
fprintf('Main stop   : %s\n', diagnostics.main_stop_reason);
fprintf('Main PG     : %.3e\n', diagnostics.main_pg_residual);
if diagnostics.entropy_enabled
    fprintf('Physical E  : %.15e\n', diagnostics.physical_energy);
    fprintf('Entropy H   : %.15e\n', diagnostics.entropy_value);
    fprintf('Augmented E : %.15e\n', diagnostics.augmented_energy);
else
    fprintf('Regularized density energy E_{epsilon,sigma}\n');
    fprintf('Fisher s    : %s\n', diagnostics.fisher_label);
    fprintf('epsilon     : %.3e\n', diagnostics.fisher_epsilon);
    fprintf('Trapping V  : %s\n', diagnostics.trapping_potential_label);
    fprintf('sigma       : %.3e\n', diagnostics.potential_sigma);
    fprintf('Potential p : %s\n', diagnostics.potential_label);
    fprintf('Target E_{eps,sigma}: %.15e\n', diagnostics.target_energy);
    fprintf('Baseline E_{eps,0}  : %.15e\n', diagnostics.baseline_energy);
end
fprintf('Window dE   : %.3e\n', diagnostics.main_energy_window_span);
fprintf('Polish      : %s\n', strrep(diagnostics.polish_solver, '_', ' '));
fprintf('Precond     : %s\n', ...
    strrep(diagnostics.polish_preconditioner, '_', ' '));
fprintf('Polish iter : %d\n', diagnostics.polish_iterations);
if diagnostics.total_pcg_z_iterations > 0 ...
        || diagnostics.total_pcg_w_iterations > 0
    fprintf('Total PCG z : %d\n', diagnostics.total_pcg_z_iterations);
    fprintf('Total PCG w : %d\n', diagnostics.total_pcg_w_iterations);
    fprintf('Max PCG iter: %d\n', diagnostics.max_pcg_iterations);
end
fprintf('PG final    : %.3e\n', diagnostics.final_pg_residual);
if isfinite(diagnostics.polish_interior_residual)
    fprintf('Interior res: %.3e\n', diagnostics.polish_interior_residual);
end
fprintf('KKT final   : %.3e\n', diagnostics.final_kkt_residual);
fprintf('dE polish   : %.3e\n', diagnostics.polish_energy_change);
fprintf('dState      : %.3e\n', diagnostics.polish_state_change);
fprintf('Mass error  : %.3e\n', diagnostics.mass_error);
fprintf('Min density : %.3e\n', diagnostics.min_density);
fprintf('Exact zeros : %d\n', diagnostics.exact_zero_count);
if isfinite(diagnostics.projected_active_count_handoff)
    fprintf('ProjAct handoff: %d (diagnostic only)\n', ...
        diagnostics.projected_active_count_handoff);
end
if isfinite(diagnostics.projected_active_count_final)
    fprintf('ProjAct final  : %d (diagnostic only)\n', ...
        diagnostics.projected_active_count_final);
end
end
