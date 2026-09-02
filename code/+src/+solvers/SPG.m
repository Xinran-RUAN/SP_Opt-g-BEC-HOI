function [rho, diagnostics, history] = SPG(problem, rho0, solver)
%SPG Spectral projected gradient with BB steps and GLL line search.
%
% Here "spectral" refers to the Barzilai-Borwein step length, not to the
% Fourier pseudospectral spatial discretization.

startTime = tic;
rho = rho0(:);
energy = src.discretization.ps.Energy(rho, problem);
gradient = src.discretization.ps.Gradient(rho, problem);
state = src.solvers.EvaluateState(rho, problem, solver, energy, gradient);
history = allocateHistory(solver.max_iter);
energyMemory = energy;
rhoPrevious = [];
gradientPrevious = [];
alpha = clampAlpha(solver.spg_alpha_reset, solver);
relativeEnergyChange = NaN;
converged = state.pg_residual <= solver.pg_tol;
failed = false;
failureMessage = '';
iterationCount = 0;

for iteration = 1:solver.max_iter
    if converged
        break;
    end
    if ~isempty(rhoPrevious)
        s = rho - rhoPrevious;
        y = gradient - gradientPrevious;
        alpha = spectralStep(s, y, iteration, problem.grid.h, solver);
    end

    [direction, descent] = projectedDirection( ...
        rho, gradient, alpha, problem, solver);
    firstMemoryIndex = max(1, numel(energyMemory) ...
        - solver.spg_nonmonotone_M + 1);
    referenceEnergy = max(energyMemory(firstMemoryIndex:end));
    [accepted, candidate, candidateEnergy, lambda, backtracks] = ...
        lineSearch(rho, direction, energy, descent, ...
        referenceEnergy, problem, solver);

    if ~accepted
        alpha = clampAlpha(solver.spg_alpha_reset, solver);
        [direction, descent] = projectedDirection( ...
            rho, gradient, alpha, problem, solver);
        [accepted, candidate, candidateEnergy, lambda, backtracksReset] = ...
            lineSearch(rho, direction, energy, descent, ...
            referenceEnergy, problem, solver);
        backtracks = backtracks + backtracksReset;
    end
    if ~accepted
        failed = true;
        failureMessage = ['SPG nonmonotone line search failed after a BB ' ...
            'reset and conservative projected-gradient retry.'];
        break;
    end

    massError = abs(src.constraints.Mass(candidate, problem.grid.h) - problem.mass);
    if min(candidate) < -1e-13 || massError > 1e-12
        error('src:solvers:SPG:FeasibilityFailure', ...
            'Line-segment trial lost simplex feasibility.');
    end
    candidateGradient = src.discretization.ps.Gradient(candidate, problem);
    candidateState = src.solvers.EvaluateState( ...
        candidate, problem, solver, candidateEnergy, candidateGradient);
    relativeEnergyChange = abs(candidateEnergy - energy) / max(1, abs(energy));

    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.energy(iteration) = candidateEnergy;
    history.pg_residual(iteration) = candidateState.pg_residual;
    history.relative_energy_change(iteration) = relativeEnergyChange;
    history.mass_error(iteration) = candidateState.mass_error;
    history.min_density(iteration) = candidateState.min_density;
    history.alpha_bb(iteration) = alpha;
    history.lambda(iteration) = lambda;
    history.backtracks(iteration) = backtracks;
    history.nonmonotone_reference(iteration) = referenceEnergy;
    history.descent(iteration) = descent;

    rhoPrevious = rho;
    gradientPrevious = gradient;
    rho = candidate;
    energy = candidateEnergy;
    gradient = candidateGradient;
    state = candidateState;
    energyMemory(end+1, 1) = energy; %#ok<AGROW>

    if solver.display && (iteration == 1 || mod(iteration, solver.display_interval) == 0)
        fprintf('%7d  %.15e  %.3e  %.3e  %.3e  %.3e  %3d\n', ...
            iteration, energy, state.pg_residual, relativeEnergyChange, ...
            alpha, lambda, backtracks);
    end
    if state.pg_residual <= solver.pg_tol
        converged = true;
    end
end

history = trimHistory(history, iterationCount);
diagnostics.iterations = iterationCount;
diagnostics.energy = energy;
diagnostics.pg_residual = state.pg_residual;
diagnostics.kkt_residual = state.kkt_residual;
diagnostics.relative_energy_change = relativeEnergyChange;
diagnostics.mass_error = state.mass_error;
diagnostics.min_density = state.min_density;
diagnostics.alpha_bb = alpha;
diagnostics.converged = converged;
diagnostics.failed = failed;
diagnostics.failure_message = failureMessage;
diagnostics.elapsed_time = toc(startTime);
diagnostics.energy_plateau = relativeEnergyChange <= solver.energy_tol;
end

function alpha = spectralStep(s, y, iteration, h, solver)
ss = h * sum(s .* s);
sy = h * sum(s .* y);
yy = h * sum(y .* y);
if sy <= solver.spg_curvature_tol || ~isfinite(sy)
    alpha = clampAlpha(solver.spg_alpha_reset, solver);
    return;
end
switch lower(char(solver.spg_bb_type))
    case 'bb1'
        alpha = ss / sy;
    case 'bb2'
        alpha = sy / yy;
    case 'alternate'
        if mod(iteration, 2) == 0
            alpha = ss / sy;
        else
            alpha = sy / yy;
        end
    otherwise
        error('src:solvers:SPG:UnknownBBType', ...
            'spg_bb_type must be bb1, bb2, or alternate.');
end
if ~isfinite(alpha) || alpha <= 0
    alpha = clampAlpha(solver.spg_alpha_reset, solver);
end
alpha = clampAlpha(alpha, solver);
end

function alpha = clampAlpha(alpha, solver)
alpha = min(solver.spg_alpha_max, max(solver.spg_alpha_min, alpha));
end

function [direction, descent] = projectedDirection( ...
    rho, gradient, alpha, problem, solver)
projected = src.solvers.ProxStep( ...
    rho, gradient, 1 / alpha, problem, solver);
direction = projected - rho;
descent = problem.grid.h * sum(gradient .* direction);
gradientNorm = sqrt(problem.grid.h * sum(gradient .^ 2));
directionNorm = sqrt(problem.grid.h * sum(direction .^ 2));
allowedPositive = solver.spg_descent_tol * max(1, gradientNorm * directionNorm);
if descent > allowedPositive
    error('src:solvers:SPG:NonDescentDirection', ...
        'Projected direction has positive h-inner-product %.3e.', descent);
end
end

function [accepted, candidate, candidateEnergy, lambda, backtracks] = ...
    lineSearch(rho, direction, currentEnergy, descent, ...
    referenceEnergy, problem, solver)
lambda = 1;
accepted = false;
candidate = rho;
candidateEnergy = currentEnergy;
backtracks = 0;
for trial = 0:solver.max_backtracks
    if lambda < solver.spg_min_lambda
        return;
    end
    candidate = rho + lambda * direction;
    candidateEnergy = src.discretization.ps.Energy(candidate, problem);
    armijoBound = referenceEnergy + solver.spg_c1 * lambda * descent;
    roundoff = 100 * eps(max([1, abs(referenceEnergy), abs(candidateEnergy)]));
    if candidateEnergy <= armijoBound + roundoff
        accepted = true;
        return;
    end
    lambda = solver.spg_backtrack * lambda;
    backtracks = backtracks + 1;
end
end

function history = allocateHistory(maxIter)
fields = {'iteration', 'elapsed_time', 'energy', 'pg_residual', ...
    'relative_energy_change', 'mass_error', 'min_density', 'alpha_bb', ...
    'lambda', 'backtracks', 'nonmonotone_reference', 'descent'};
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
