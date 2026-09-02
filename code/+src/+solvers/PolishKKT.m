function result = PolishKKT(rho0, problem, options)
%POLISHKKT Dispatch complete-target polishing by the detected active set.

options = applyDispatcherDefaults(options, problem.grid.N);
initialGradient = src.discretization.ps.Gradient(rho0, problem);
initialDiagnostic = src.solvers.InteriorKKTResidual( ...
    rho0, problem, options, initialGradient);
requestedSolver = char(options.linear_solver);

if strcmpi(requestedSolver, 'pdas_gmres')
    result = runPDASGMRES(rho0, problem, pdasOptions(options));
    result.polish_solver = 'pdas_gmres';
    result.switched_to_pdas = false;
    result.fallback_reason = '';
    result = addPcgDefaults(result, options.preconditioner);
    result = attachDispatchDiagnostics(result, initialDiagnostic);
    return;
end

[interiorAllowed, eligibility] = interiorEligibility(rho0, problem, options);
if interiorAllowed
    if options.display
        fprintf(['Interior Newton-PCG polish\n' ...
            'iter   int_res    full_PG    PCG(z/w)   ' ...
            'relres(z/w)       step   projAct\n']);
    end
    interior = src.solvers.PolishInteriorNewtonPCG( ...
        rho0, problem, options);
    if ~interior.request_pdas_fallback
        result = attachDispatchDiagnostics(interior, initialDiagnostic);
        return;
    end
    boundaryEvidence = any(interior.rho == 0) ...
        || ~eligibility.vacuum_compatible;
    if ~boundaryEvidence
        interior.failed = true;
        interior.failure_message = sprintf([ ...
            '%s Projected active count %d is diagnostic only; ' ...
            'PDAS fallback suppressed without exact boundary evidence.'], ...
            interior.fallback_reason, interior.projected_active_count);
        interior.status = 'interior_failure_without_boundary_evidence';
        interior.request_pdas_fallback = false;
        result = attachDispatchDiagnostics(interior, initialDiagnostic);
        return;
    end
    fallbackReason = interior.fallback_reason;
    fallbackState = interior.rho;
else
    fallbackReason = eligibility.reason;
    fallbackState = rho0;
    interior = [];
end

if ~options.allow_pdas_fallback
    if isempty(interior)
        result = rejectedInteriorResult( ...
            rho0, problem, initialDiagnostic, fallbackReason, options);
    else
        result = interior;
        result.failed = true;
        result.failure_message = fallbackReason;
        result.status = 'pdas_fallback_disabled';
    end
    result = attachDispatchDiagnostics(result, initialDiagnostic);
    return;
end
if options.display
    fprintf('\n[Interior PCG -> PDAS-GMRES]\nreason = %s\n', ...
        fallbackReason);
end
pdas = runPDASGMRES(fallbackState, problem, pdasOptions(options));
initialEnergy = src.discretization.ps.Energy(rho0, problem);
pdas.energy_change = pdas.energy - initialEnergy;
pdas.state_change = sqrt(problem.grid.h ...
    * sum((pdas.rho - rho0(:)) .^ 2));
pdas.pg_residual_before = initialDiagnostic.full_pg_residual;
pdas.polish_solver = 'pdas_gmres_fallback';
pdas.switched_to_pdas = true;
pdas.fallback_reason = fallbackReason;
if isempty(interior)
    pdas = addPcgDefaults(pdas, options.preconditioner);
else
    combinedHistory.interior = interior.history;
    combinedHistory.pdas = pdas.history;
    pdas.history = combinedHistory;
    pdas.iterations = interior.iterations + pdas.iterations;
    pdas.elapsed_time = interior.elapsed_time + pdas.elapsed_time;
    pdas.total_pcg_z_iterations = interior.total_pcg_z_iterations;
    pdas.total_pcg_w_iterations = interior.total_pcg_w_iterations;
    pdas.max_pcg_iterations = interior.max_pcg_iterations;
    pdas.preconditioner_shift = interior.preconditioner_shift;
    pdas.preconditioner = interior.preconditioner;
end
result = pdas;
result = attachDispatchDiagnostics(result, initialDiagnostic);
end

function result = attachDispatchDiagnostics(result, diagnostic)
result.handoff_projected_active_count = ...
    diagnostic.projected_active_count;
result.handoff_projected_free_count = diagnostic.projected_free_count;
result.handoff_exact_zero_count = diagnostic.exact_zero_count;
result.handoff_min_density = diagnostic.min_density;
end

function result = runPDASGMRES(rho0, problem, options)
% Existing matrix-free PDAS-GMRES fallback for inequality-active states.

options = applyDefaults(options);
startTime = tic;
rho = rho0(:);
state = fullState(rho, problem, options);
initialState = state;
history = allocateHistory(options.max_iter);
converged = state.pg_residual <= options.pg_tol;
attempted = false;
failed = false;
failureMessage = '';
iterationCount = 0;
previousFree = [];

if state.pg_residual > options.entry_pg_tol
    failureMessage = sprintf( ...
        'KKT entry rejected: full PG %.3e exceeds entry tolerance %.3e.', ...
        state.pg_residual, options.entry_pg_tol);
    result = packageResult();
    return;
end

attempted = true;
for iteration = 1:options.max_iter
    if converged
        break;
    end
    projected = src.constraints.ProjectPositiveConservative( ...
        rho - options.active_step * state.gradient, ...
        problem.mass, problem.grid.h, options.projection_tol);
    free = projected > options.active_tol;
    if ~any(free)
        [~, freeIndex] = max(projected);
        free(freeIndex) = true;
    end
    active = ~free;
    if isempty(previousFree)
        activeSetChanges = 0;
    else
        activeSetChanges = nnz(xor(free, previousFree));
    end
    previousFree = free;

    lambda = -mean(state.gradient(free));
    activeCorrection = zeros(problem.grid.N, 1);
    activeCorrection(active) = -rho(active);
    hessianActive = src.solvers.FullHessianAction( ...
        rho, activeCorrection, problem);
    massResidual = src.constraints.Mass(rho, problem.grid.h) - problem.mass;
    rightHandSide = [ ...
        -(state.gradient(free) + lambda + hessianActive(free)); ...
        (-massResidual - problem.grid.h ...
        * sum(activeCorrection(active))) / problem.grid.h];

    operator = @(z) kktAction(z, free, rho, problem);
    numberOfUnknowns = nnz(free) + 1;
    restart = min(numberOfUnknowns, options.gmres_maxit);
    maxOuter = min(options.gmres_maxit, numberOfUnknowns);
    diagonal = fullHessianDiagonal(rho, problem);
    inverseScaling = @(v) kktInverseScaling( ...
        v, free, diagonal);
    scaledOperator = @(v) operator(inverseScaling(v));
    [scaledZ, gmresFlag, gmresRelativeResidual, gmresIteration] = gmres( ...
        scaledOperator, rightHandSide, restart, ...
        options.gmres_tol, maxOuter);
    z = inverseScaling(scaledZ);
    gmresIterations = countGmresIterations(gmresIteration, restart);
    if gmresFlag ~= 0 && gmresRelativeResidual > options.gmres_failure_tol
        failed = true;
        failureMessage = sprintf( ...
            'KKT GMRES failed: flag=%d, relative residual=%.3e.', ...
            gmresFlag, gmresRelativeResidual);
        break;
    end

    correction = activeCorrection;
    correction(free) = z(1:end-1);
    negativeDirection = correction < 0;
    if any(negativeDirection)
        maximumPositiveStep = min( ...
            -rho(negativeDirection) ./ correction(negativeDirection));
    else
        maximumPositiveStep = inf;
    end
    step = min(1, 0.99 * maximumPositiveStep);
    accepted = false;
    backtracks = 0;
    trialState = state;
    trialRho = rho;

    while step >= options.min_step && backtracks <= options.max_backtracks
        trialRho = rho + step * correction;
        [trialRho, correctionApplied] = correctMassRoundoff( ...
            trialRho, free, problem, options);
        trialMassError = abs(src.constraints.Mass( ...
            trialRho, problem.grid.h) - problem.mass);
        if min(trialRho) >= -options.roundoff_tol ...
                && trialMassError <= options.mass_tol ...
                && correctionApplied
            trialState = fullState(trialRho, problem, options);
            residualBound = (1 - options.residual_c1 * step) ...
                * state.pg_residual;
            if trialState.pg_residual <= residualBound ...
                    || trialState.pg_residual <= options.pg_tol
                accepted = true;
                break;
            end
        end
        step = options.backtrack * step;
        backtracks = backtracks + 1;
    end

    if ~accepted
        failed = true;
        failureMessage = sprintf([ ...
            'KKT residual globalization failed: PG %.3e -> %.3e, ' ...
            'initial step %.3e, final step %.3e, active/free %d/%d, ' ...
            'GMRES flag %d relres %.3e.'], ...
            state.pg_residual, trialState.pg_residual, ...
            min(1, 0.99 * maximumPositiveStep), step, ...
            nnz(active), nnz(free), gmresFlag, gmresRelativeResidual);
        break;
    end

    rho = trialRho;
    state = trialState;
    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.energy(iteration) = state.energy;
    history.pg_residual(iteration) = state.pg_residual;
    history.kkt_residual(iteration) = state.kkt_residual;
    history.mass_error(iteration) = state.mass_error;
    history.min_density(iteration) = state.min_density;
    history.active_count(iteration) = nnz(active);
    history.free_count(iteration) = nnz(free);
    history.active_set_changes(iteration) = activeSetChanges;
    history.gmres_iterations(iteration) = gmresIterations;
    history.gmres_flag(iteration) = gmresFlag;
    history.gmres_relative_residual(iteration) = gmresRelativeResidual;
    history.step(iteration) = step;
    history.backtracks(iteration) = backtracks;

    if options.display
        fprintf('polish %3d  %.3e  %.3e  %5d  %4d  %.3e\n', ...
            iteration, state.pg_residual, state.kkt_residual, ...
            nnz(active), gmresIterations, step);
    end
    converged = state.pg_residual <= options.pg_tol;
end

result = packageResult();

    function output = packageResult()
        output.rho = rho;
        output.energy = state.energy;
        output.history = trimHistory(history, iterationCount);
        output.polish_attempted = attempted;
        output.polish_converged = converged;
        output.failed = failed;
        output.failure_message = failureMessage;
        if converged
            output.status = 'converged';
        elseif failed
            output.status = 'failed';
        elseif attempted
            output.status = 'max_iterations';
        else
            output.status = 'entry_rejected';
        end
        output.iterations = iterationCount;
        output.elapsed_time = toc(startTime);
        output.pg_residual_before = initialState.pg_residual;
        output.pg_residual = state.pg_residual;
        output.kkt_residual = state.kkt_residual;
        output.mass_error = state.mass_error;
        output.min_density = state.min_density;
        output.active_count = state.active_count;
        output.free_count = state.free_count;
        output.energy_change = state.energy - initialState.energy;
        output.state_change = sqrt(problem.grid.h ...
            * sum((rho - rho0(:)) .^ 2));
        output.polish_solver = 'pdas_gmres';
        output.preconditioner = 'kkt_diagonal_scaling';
        output.total_gmres_iterations = ...
            sum(output.history.gmres_iterations);
        if isempty(output.history.gmres_iterations)
            output.max_gmres_iterations = 0;
        else
            output.max_gmres_iterations = ...
                max(output.history.gmres_iterations);
        end
        output.request_pdas_fallback = false;
        output.switched_to_pdas = false;
        output.fallback_reason = '';
    end
end

function state = fullState(rho, problem, options)
state.energy = src.discretization.ps.Energy(rho, problem);
state.gradient = src.discretization.ps.Gradient(rho, problem);
[state.pg_residual, state.pg_mapping, state.projected_state] = ...
    src.solvers.FullGradientMapping(rho, state.gradient, problem, options);
state.kkt = src.solvers.KKTResidual( ...
    rho, state.gradient, problem, options.active_tol);
state.kkt_residual = state.kkt.kkt_residual;
state.mass_error = abs(src.constraints.Mass( ...
    rho, problem.grid.h) - problem.mass);
state.min_density = min(rho);
state.active_count = state.kkt.active_count;
state.free_count = state.kkt.free_count;
end

function value = kktAction(z, free, rho, problem)
numberOfFree = nnz(free);
direction = zeros(size(rho));
direction(free) = z(1:numberOfFree);
hessianDirection = src.solvers.FullHessianAction(rho, direction, problem);
value = [hessianDirection(free) + z(end); ...
    sum(direction(free))];
end

function diagonal = fullHessianDiagonal(rho, problem)
% Exact nodal diagonal assembled without forming the dense Fourier Hessian.
basis = zeros(problem.grid.N, 1);
basis(1) = 1;
derivativeColumn = src.discretization.ps.FirstDerivative( ...
    basis, problem.plan);
[s, ds, d2s] = src.regularization.EvaluateFisher( ...
    rho, src.regularization.ResolveFisher(problem));
q = src.discretization.ps.FirstDerivative(rho, problem.plan);
weight = 1 ./ s;
fisherDifferentialDiagonal = 0.25 * real(ifft( ...
    conj(fft(derivativeColumn .^ 2)) .* fft(weight)));
denominatorDerivative = d2s ./ (s .^ 2) ...
    - 2 * (ds .^ 2) ./ (s .^ 3);
fisherLocalDiagonal = -0.125 * q .^ 2 .* denominatorDerivative;
dtDDiagonal = sum(derivativeColumn .^ 2);
diagonal = fisherDifferentialDiagonal + fisherLocalDiagonal ...
    + problem.beta + problem.delta * dtDDiagonal;
potential = problem.potential_regularization;
fisher = src.regularization.ResolveFisher(problem);
[~, ~, potentialSecondDerivative, ~] = src.potential.Evaluate( ...
    rho, potential, fisher.epsilon);
diagonal = diagonal + problem.V(:) .* potentialSecondDerivative;
diagonal = max(diagonal, problem.beta);
end

function value = kktInverseScaling(z, free, diagonal)
freeDiagonal = diagonal(free);
schurScale = max(sum(1 ./ freeDiagonal), eps);
value = [z(1:end-1) ./ freeDiagonal; z(end) / schurScale];
end

function count = countGmresIterations(iteration, restart)
if iteration(1) == 0
    count = iteration(2);
else
    count = (iteration(1) - 1) * restart + iteration(2);
end
end

function [rho, accepted] = correctMassRoundoff(rho, free, problem, options)
massResidual = src.constraints.Mass(rho, problem.grid.h) - problem.mass;
if abs(massResidual) <= options.mass_tol
    accepted = true;
    if massResidual ~= 0
        rho(free) = rho(free) ...
            - massResidual / (problem.grid.h * nnz(free));
    end
else
    accepted = false;
end
end

function options = applyDefaults(options)
if nargin < 1 || isempty(options)
    options = struct();
end
if isfield(options, 'tau_active') && ~isfield(options, 'active_step')
    options.active_step = options.tau_active;
end
defaults.entry_pg_tol = 1e-5;
defaults.pg_tol = 1e-12;
defaults.vacuum_slope_tol = 1e-14;
defaults.convexity_tol = 1e-14;
defaults.max_iter = 50;
defaults.active_step = 1;
defaults.active_tol = 1e-12;
defaults.gmres_tol = 1e-10;
defaults.gmres_maxit = 200;
defaults.gmres_failure_tol = 0.5;
defaults.residual_step = 1;
defaults.projection_tol = 1e-14;
defaults.residual_c1 = 1e-4;
defaults.backtrack = 0.5;
defaults.min_step = 1e-12;
defaults.max_backtracks = 30;
defaults.mass_tol = 1e-12;
defaults.roundoff_tol = 1e-13;
defaults.display = true;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
end

function history = allocateHistory(maxIter)
fields = {'iteration', 'elapsed_time', 'energy', 'pg_residual', ...
    'kkt_residual', 'mass_error', 'min_density', 'active_count', ...
    'free_count', 'active_set_changes', 'gmres_iterations', 'gmres_flag', ...
    'gmres_relative_residual', 'step', 'backtracks'};
for j = 1:numel(fields)
    history.(fields{j}) = nan(maxIter, 1);
end
end

function history = trimHistory(history, count)
fields = fieldnames(history);
for j = 1:numel(fields)
    history.(fields{j}) = history.(fields{j})(1:count, :);
end
end

function options = applyDispatcherDefaults(options, N)
if nargin < 1 || isempty(options)
    options = struct();
end
defaults.linear_solver = 'interior_pcg_schur';
defaults.preconditioner = 'fd_variable';
defaults.allow_pdas_fallback = true;
defaults.pdas_max_iter = 50;
defaults.max_iter = 20;
defaults.pcg_tol_max = 1e-2;
defaults.pcg_tol_min = 1e-10;
defaults.pcg_forcing_factor = 0.1;
defaults.pcg_maxit = min(500, max(100, N));
defaults.fraction_to_boundary = 0.995;
defaults.residual_armijo = 1e-4;
defaults.stagnation_window = 5;
defaults.stagnation_rel_improvement = 1e-2;
defaults.acceptable_floor = 1e-10;
defaults.active_step = 1;
defaults.active_tol = 1e-12;
defaults.projection_tol = 1e-14;
defaults.residual_step = 1;
defaults.pg_tol = 1e-12;
defaults.vacuum_slope_tol = 1e-14;
defaults.convexity_tol = 1e-14;
defaults.display = true;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
if ~ismember(lower(char(options.linear_solver)), ...
        {'interior_pcg_schur', 'auto', 'pdas_gmres'})
    error('src:solvers:PolishKKT:LinearSolver', ...
        'linear_solver must be interior_pcg_schur, auto, or pdas_gmres.');
end
end

function tf = knownSpdRegime(rho, problem)
fisher = src.regularization.ResolveFisher(problem);
[r, ~, d2r] = src.regularization.EvaluateFisher(rho, fisher);
potential = problem.potential_regularization;
d2p = potential.d2p_sigma(rho);
tf = problem.beta > 0 && problem.delta >= 0 ...
    && min(problem.V(:)) >= -1e-14 && all(r > 0) ...
    && all(d2r <= 1e-13) && all(d2p >= -1e-13);
end

function [eligible, diagnostic] = interiorEligibility(rho, problem, options)
fisher = src.regularization.ResolveFisher(problem);
r0 = fisher.r_epsilon(zeros(size(rho(1))));
dp0 = problem.potential_regularization.dp_sigma(0);
diagnostic.strictly_positive = all(isfinite(rho)) && all(rho > 0);
diagnostic.vacuum_compatible = isscalar(r0) && isfinite(r0) && r0 > 0 ...
    && isscalar(dp0) && isfinite(dp0) ...
    && abs(dp0) <= options.vacuum_slope_tol;
diagnostic.strong_convex_regime = problem.beta > 0 ...
    && problem.delta >= 0 ...
    && min(problem.V(:)) >= -options.convexity_tol ...
    && knownSpdRegime(rho, problem);
eligible = diagnostic.strictly_positive ...
    && diagnostic.vacuum_compatible ...
    && diagnostic.strong_convex_regime;
if ~diagnostic.strictly_positive
    diagnostic.reason = 'exact zero/nonpositive density requires boundary handling';
elseif ~diagnostic.vacuum_compatible
    diagnostic.reason = 'potential is not vacuum-compatible for interior polish';
else
    diagnostic.reason = 'known strong-convexity/SPD guard is not satisfied';
end
end

function options = pdasOptions(options)
options.max_iter = options.pdas_max_iter;
end

function result = addPcgDefaults(result, preconditioner)
result.total_pcg_z_iterations = 0;
result.total_pcg_w_iterations = 0;
result.max_pcg_iterations = 0;
result.preconditioner_shift = 0;
if ~isfield(result, 'preconditioner') || isempty(result.preconditioner)
    result.preconditioner = preconditioner;
end
if ~isfield(result, 'interior_residual')
    result.interior_residual = NaN;
end
end

function result = rejectedInteriorResult( ...
    rho, problem, diagnostic, reason, options)
result.rho = rho(:);
result.energy = src.discretization.ps.Energy(rho, problem);
result.history = struct();
result.polish_attempted = false;
result.polish_converged = false;
result.failed = true;
result.failure_message = reason;
result.status = 'pdas_fallback_disabled';
result.iterations = 0;
result.elapsed_time = 0;
result.pg_residual_before = diagnostic.full_pg_residual;
result.pg_residual = diagnostic.full_pg_residual;
result.full_pg_residual = diagnostic.full_pg_residual;
result.interior_residual = diagnostic.interior_residual;
kkt = src.solvers.KKTResidual( ...
    rho, diagnostic.gradient, problem, options.active_tol);
result.kkt_residual = kkt.kkt_residual;
result.mass_error = diagnostic.mass_error;
result.min_density = diagnostic.min_density;
result.active_count = diagnostic.active_count;
result.free_count = diagnostic.free_count;
result.energy_change = 0;
result.state_change = 0;
result.polish_solver = 'none';
result.switched_to_pdas = false;
result.fallback_reason = reason;
result.request_pdas_fallback = false;
result = addPcgDefaults(result, options.preconditioner);
end
