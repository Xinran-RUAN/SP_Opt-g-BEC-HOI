function result = PolishInteriorNewtonPCG(rho0, problem, options)
%POLISHINTERIORNEWTONPCG Equality Newton using SPD PCG-Schur solves.
%
% The interior target Hessian is SPD in the current strongly convex regime.
% PCG is therefore used instead of GMRES on the indefinite KKT saddle
% system. All Hessian products still use the spectral matrix-free action.

options = applyDefaults(options, problem.grid.N);
startTime = tic;
rho = rho0(:);
energy = src.discretization.ps.Energy(rho, problem);
state = src.solvers.InteriorKKTResidual(rho, problem, options);
initialEnergy = energy;
initialState = state;
history = allocateHistory(options.max_iter);
iterationCount = 0;
converged = state.full_pg_residual <= options.pg_tol;
failed = false;
failureMessage = '';
status = 'not_attempted';
requestPdasFallback = false;
fallbackReason = '';
firstDirection = [];
firstDeltaLambda = NaN;
totalPcgZ = 0;
totalPcgW = 0;
maximumPcg = 0;
maximumShift = 0;

for iteration = 1:options.max_iter
    if converged
        status = 'converged';
        break;
    end
    shift = 0;
    Hfun = @(v) src.solvers.FullHessianAction(rho, v, problem);
    if strcmpi(options.preconditioner, 'fd_variable')
        [P, preconditionerInfo] = ...
            src.solvers.BuildFDHessianPreconditioner(rho, problem);
        [R, cholFlag, shift] = factorPreconditioner(P, preconditionerInfo);
        maximumShift = max(maximumShift, shift);
        if cholFlag ~= 0
            requestPdasFallback = true;
            fallbackReason = sprintf([ ...
                'FD preconditioner Cholesky failed: min diagonal %.3e, ' ...
                'symmetry error %.3e, machine shift %.3e.'], ...
                preconditionerInfo.min_diagonal, ...
                preconditionerInfo.symmetry_error, shift);
            status = 'pcg_failure';
            break;
        end
        Mfun = @(v) R' \ (R \ v);
    else
        Mfun = [];
    end
    pcgTolerance = min(options.pcg_tol_max, ...
        max(options.pcg_tol_min, options.pcg_forcing_factor ...
        * sqrt(state.interior_residual)));
    [z, zFlag, zRelres, zIterations] = solvePcg( ...
        Hfun, state.stationarity, pcgTolerance, options.pcg_maxit, Mfun);
    [w, wFlag, wRelres, wIterations] = solvePcg( ...
        Hfun, ones(problem.grid.N, 1), pcgTolerance, ...
        options.pcg_maxit, Mfun);
    if zFlag ~= 0 || wFlag ~= 0
        retryTolerance = min(10 * pcgTolerance, 1e-4);
        if zFlag ~= 0
            [z, zFlag, zRelres, zIterations] = solvePcg( ...
                Hfun, state.stationarity, retryTolerance, ...
                options.pcg_maxit, Mfun);
        end
        if wFlag ~= 0
            [w, wFlag, wRelres, wIterations] = solvePcg( ...
                Hfun, ones(problem.grid.N, 1), retryTolerance, ...
                options.pcg_maxit, Mfun);
        end
    end
    totalPcgZ = totalPcgZ + zIterations;
    totalPcgW = totalPcgW + wIterations;
    maximumPcg = max([maximumPcg, zIterations, wIterations]);
    if zFlag ~= 0 || wFlag ~= 0
        requestPdasFallback = true;
        fallbackReason = sprintf([ ...
            'PCG failure: z flag/relres %d/%.3e, ' ...
            'w flag/relres %d/%.3e.'], ...
            zFlag, zRelres, wFlag, wRelres);
        status = 'pcg_failure';
        break;
    end

    h = problem.grid.h;
    denominator = h * sum(w);
    if ~isfinite(denominator) || denominator <= 0
        requestPdasFallback = true;
        fallbackReason = 'Nonpositive PCG-Schur denominator.';
        status = 'pcg_failure';
        break;
    end
    deltaLambda = (state.mass_residual - h * sum(z)) / denominator;
    direction = -z - deltaLambda * w;
    massLinearizationResidual = abs( ...
        h * sum(direction) + state.mass_residual);
    if iteration == 1
        firstDirection = direction;
        firstDeltaLambda = deltaLambda;
    end

    negative = direction < 0;
    if any(negative)
        maximumStep = min(-rho(negative) ./ direction(negative));
        step = min(1, options.fraction_to_boundary * maximumStep);
    else
        step = 1;
    end
    accepted = false;
    backtracks = 0;
    trialRho = rho;
    trialState = state;
    trialEnergy = energy;
    while step >= options.min_step && backtracks <= options.max_backtracks
        trialRho = rho + step * direction;
        [trialRho, massCorrectionOk] = correctMassRoundoff( ...
            trialRho, problem, options);
        if massCorrectionOk && min(trialRho) > 0
            trialGradient = src.discretization.ps.Gradient(trialRho, problem);
            trialState = src.solvers.InteriorKKTResidual( ...
                trialRho, problem, options, trialGradient);
            if trialState.interior_residual <= ...
                    (1 - options.residual_armijo * step) ...
                    * state.interior_residual ...
                    || trialState.full_pg_residual <= options.pg_tol
                trialEnergy = src.discretization.ps.Energy(trialRho, problem);
                accepted = true;
                break;
            end
        end
        step = options.backtrack * step;
        backtracks = backtracks + 1;
    end
    if ~accepted
        requestPdasFallback = true;
        fallbackReason = 'Interior residual line search failed.';
        status = 'residual_line_search_failure';
        break;
    end

    rho = trialRho;
    state = trialState;
    energy = trialEnergy;
    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.energy(iteration) = energy;
    history.interior_residual(iteration) = state.interior_residual;
    history.full_pg_residual(iteration) = state.full_pg_residual;
    history.mass_error(iteration) = state.mass_error;
    history.min_density(iteration) = state.min_density;
    history.active_count(iteration) = state.projected_active_count;
    history.projected_active_count(iteration) = ...
        state.projected_active_count;
    history.exact_zero_count(iteration) = state.exact_zero_count;
    history.pcg_tolerance(iteration) = pcgTolerance;
    history.pcg_z_flag(iteration) = zFlag;
    history.pcg_z_relres(iteration) = zRelres;
    history.pcg_z_iter(iteration) = zIterations;
    history.pcg_w_flag(iteration) = wFlag;
    history.pcg_w_relres(iteration) = wRelres;
    history.pcg_w_iter(iteration) = wIterations;
    history.mass_linearization_residual(iteration) = ...
        massLinearizationResidual;
    history.step(iteration) = step;
    history.backtracks(iteration) = backtracks;
    history.preconditioner_shift(iteration) = shift;

    if options.display
        fprintf(['polish %2d  %.3e  %.3e  %4d/%-4d  ' ...
            '%.1e/%.1e  %.3e  %d\n'], iteration, ...
            state.interior_residual, state.full_pg_residual, ...
            zIterations, wIterations, zRelres, wRelres, step, ...
            state.projected_active_count);
    end
    converged = state.full_pg_residual <= options.pg_tol;
    if converged
        status = 'converged';
        break;
    end
    if iterationCount > options.stagnation_window
        oldResidual = history.full_pg_residual( ...
            iterationCount - options.stagnation_window);
        relativeImprovement = (oldResidual - state.full_pg_residual) ...
            / max(oldResidual, realmin);
        if relativeImprovement < options.stagnation_rel_improvement
            if state.full_pg_residual <= options.acceptable_floor
                status = 'residual_floor';
            else
                status = 'stagnated_above_tolerance';
                failed = true;
                failureMessage = sprintf( ...
                    'Interior polish stagnated at full PG %.3e.', ...
                    state.full_pg_residual);
            end
            break;
        end
    end
end

if strcmp(status, 'not_attempted')
    if converged
        status = 'converged';
    else
        status = 'max_iterations';
    end
end
kkt = src.solvers.KKTResidual( ...
    rho, state.gradient, problem, options.active_tol);
result.rho = rho;
result.energy = energy;
result.history = trimHistory(history, iterationCount);
result.polish_attempted = true;
result.polish_converged = converged;
result.failed = failed;
result.failure_message = failureMessage;
result.status = status;
result.iterations = iterationCount;
result.elapsed_time = toc(startTime);
result.pg_residual_before = initialState.full_pg_residual;
result.pg_residual = state.full_pg_residual;
result.full_pg_residual = state.full_pg_residual;
result.interior_residual = state.interior_residual;
result.kkt_residual = kkt.kkt_residual;
result.mass_error = state.mass_error;
result.min_density = state.min_density;
result.active_count = state.projected_active_count;
result.free_count = state.projected_free_count;
result.projected_active_count = state.projected_active_count;
result.projected_free_count = state.projected_free_count;
result.exact_zero_count = state.exact_zero_count;
result.strict_positive_count = state.strict_positive_count;
result.energy_change = energy - initialEnergy;
result.state_change = sqrt(problem.grid.h * sum((rho - rho0(:)) .^ 2));
result.polish_solver = 'interior_newton_pcg';
result.preconditioner = options.preconditioner;
result.total_pcg_z_iterations = totalPcgZ;
result.total_pcg_w_iterations = totalPcgW;
result.max_pcg_iterations = maximumPcg;
result.total_gmres_iterations = 0;
result.max_gmres_iterations = 0;
result.preconditioner_shift = maximumShift;
result.request_pdas_fallback = requestPdasFallback;
result.fallback_reason = fallbackReason;
result.switched_to_pdas = false;
result.first_newton_direction = firstDirection;
result.first_delta_lambda = firstDeltaLambda;
end

function [R, flag, shift] = factorPreconditioner(P, info)
shift = 0;
if isfield(info, 'factorization') && strcmpi(info.factorization, 'ichol')
    try
        setup.type = 'ict';
        setup.droptol = 1e-3;
        setup.diagcomp = 0;
        R = ichol(P, setup);
        flag = 0;
    catch
        shift = 100 * eps(class(P)) * max(1, norm(P, 1));
        try
            setup.diagcomp = shift;
            R = ichol(P, setup);
            flag = 0;
        catch
            R = sparse(size(P, 1), size(P, 2));
            flag = 1;
        end
    end
else
    [R, flag] = chol(P, 'lower');
    if flag ~= 0
        shift = 100 * eps(class(P)) * max(1, norm(P, 1));
        [R, flag] = chol(P + shift * speye(size(P, 1)), 'lower');
    end
end
end

function [x, flag, relativeResidual, iterations] = ...
    solvePcg(operator, rhs, tolerance, maxIterations, preconditioner)
if isempty(preconditioner)
    [x, flag, relativeResidual, iterations] = pcg( ...
        operator, rhs, tolerance, maxIterations);
else
    [x, flag, relativeResidual, iterations] = pcg( ...
        operator, rhs, tolerance, maxIterations, preconditioner);
end
end

function [rho, accepted] = correctMassRoundoff(rho, problem, options)
massResidual = src.constraints.Mass(rho, problem.grid.h) - problem.mass;
accepted = abs(massResidual) <= options.mass_tol;
if accepted && massResidual ~= 0
    rho = rho - massResidual ...
        / (problem.grid.h * problem.grid.N);
end
end

function options = applyDefaults(options, N)
defaults.pg_tol = 1e-12;
defaults.max_iter = 20;
defaults.active_step = 1;
defaults.active_tol = 1e-12;
defaults.projection_tol = 1e-14;
defaults.residual_step = 1;
defaults.pcg_tol_max = 1e-2;
defaults.pcg_tol_min = 1e-10;
defaults.pcg_forcing_factor = 0.1;
defaults.pcg_maxit = min(500, max(100, N));
defaults.preconditioner = 'fd_variable';
defaults.fraction_to_boundary = 0.995;
defaults.residual_armijo = 1e-4;
defaults.backtrack = 0.5;
defaults.min_step = 1e-10;
defaults.max_backtracks = 40;
defaults.mass_tol = 1e-12;
defaults.roundoff_tol = 1e-13;
defaults.stagnation_window = 5;
defaults.stagnation_rel_improvement = 1e-2;
defaults.acceptable_floor = 1e-10;
defaults.display = true;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
if ~ismember(lower(char(options.preconditioner)), {'fd_variable', 'none'})
    error('src:solvers:PolishInteriorNewtonPCG:Preconditioner', ...
        'preconditioner must be fd_variable or none.');
end
end

function history = allocateHistory(maxIterations)
fields = {'iteration', 'elapsed_time', 'energy', 'interior_residual', ...
    'full_pg_residual', 'mass_error', 'min_density', 'active_count', ...
    'projected_active_count', 'exact_zero_count', ...
    'pcg_tolerance', 'pcg_z_flag', 'pcg_z_relres', 'pcg_z_iter', ...
    'pcg_w_flag', 'pcg_w_relres', 'pcg_w_iter', ...
    'mass_linearization_residual', 'step', 'backtracks', ...
    'preconditioner_shift'};
for j = 1:numel(fields)
    history.(fields{j}) = nan(maxIterations, 1);
end
end

function history = trimHistory(history, count)
fields = fieldnames(history);
for j = 1:numel(fields)
    history.(fields{j}) = history.(fields{j})(1:count, :);
end
end
