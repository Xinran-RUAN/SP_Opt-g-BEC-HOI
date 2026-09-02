function [rho, diagnostics, history] = ISTA(problem, rho0, solver)
%ISTA Reference projected-gradient solver with monotone backtracking.

startTime = tic;
rho = rho0(:);
energy = src.discretization.ps.Energy(rho, problem);
gradient = src.discretization.ps.Gradient(rho, problem);
state = src.solvers.EvaluateState(rho, problem, solver, energy, gradient);
objective = state.augmented_energy;
L = solver.L0;
history = allocateHistory(solver.max_iter);
converged = state.pg_residual <= solver.pg_tol;
relativeEnergyChange = NaN;
iterationCount = 0;

for iteration = 1:solver.max_iter
    if converged
        break;
    end
    trialL = max(solver.L0, solver.ista_L_decrease * L);
    [rhoNew, energyNew, L] = backtrackingStep( ...
        rho, energy, gradient, trialL, problem, solver);
    gradientNew = src.discretization.ps.Gradient(rhoNew, problem);
    stateNew = src.solvers.EvaluateState( ...
        rhoNew, problem, solver, energyNew, gradientNew);
    objectiveNew = stateNew.augmented_energy;
    relativeEnergyChange = abs(objectiveNew - objective) ...
        / max(1, abs(objective));

    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.energy(iteration) = energyNew;
    history.physical_energy(iteration) = stateNew.physical_energy;
    history.entropy_value(iteration) = stateNew.entropy_value;
    history.augmented_energy(iteration) = stateNew.augmented_energy;
    history.pg_residual(iteration) = stateNew.pg_residual;
    history.relative_energy_change(iteration) = relativeEnergyChange;
    history.mass_error(iteration) = stateNew.mass_error;
    history.min_density(iteration) = stateNew.min_density;
    history.accepted_L(iteration) = L;
    history.restart(iteration) = false;

    rho = rhoNew;
    energy = energyNew;
    objective = objectiveNew;
    gradient = gradientNew;
    state = stateNew;
    if solver.display && (iteration == 1 || mod(iteration, solver.display_interval) == 0)
        fprintf('%7d  %.15e  %.3e  %.3e  %.3e  %d\n', ...
            iteration, energy, state.pg_residual, relativeEnergyChange, L, 0);
    end
    if state.pg_residual <= solver.pg_tol
        converged = true;
    end
end

history = trimHistory(history, iterationCount);
diagnostics.iterations = iterationCount;
diagnostics.energy = energy;
diagnostics.physical_energy = state.physical_energy;
diagnostics.entropy_value = state.entropy_value;
diagnostics.augmented_energy = state.augmented_energy;
diagnostics.pg_residual = state.pg_residual;
diagnostics.kkt_residual = state.kkt_residual;
diagnostics.relative_energy_change = relativeEnergyChange;
diagnostics.mass_error = state.mass_error;
diagnostics.min_density = state.min_density;
diagnostics.accepted_L = L;
diagnostics.converged = converged;
diagnostics.failed = false;
diagnostics.failure_message = '';
diagnostics.restarts = 0;
diagnostics.elapsed_time = toc(startTime);
diagnostics.energy_plateau = relativeEnergyChange <= solver.energy_tol;
end

function [candidate, candidateEnergy, acceptedL] = backtrackingStep( ...
    rho, energy, gradient, initialL, problem, solver)
acceptedL = initialL;
for trial = 1:solver.max_backtracks
    candidate = src.solvers.ProxStep( ...
        rho, gradient, acceptedL, problem, solver);
    candidateEnergy = src.discretization.ps.Energy(candidate, problem);
    difference = candidate - rho;
    majorant = energy ...
        + problem.grid.h * sum(gradient .* difference) ...
        + 0.5 * acceptedL * problem.grid.h * sum(difference .^ 2);
    roundoff = 100 * eps(max([1, abs(energy), abs(candidateEnergy)]));
    if candidateEnergy <= majorant + roundoff
        return;
    end
    acceptedL = solver.backtrack_factor * acceptedL;
    if ~isfinite(acceptedL)
        break;
    end
end
error('src:solvers:ISTA:BacktrackingFailure', ...
    'Backtracking failed to establish the h-inner-product majorization.');
end

function history = allocateHistory(maxIter)
fields = {'iteration', 'elapsed_time', 'energy', 'physical_energy', ...
    'entropy_value', 'augmented_energy', 'pg_residual', ...
    'relative_energy_change', 'mass_error', 'min_density', 'accepted_L'};
for j = 1:numel(fields)
    history.(fields{j}) = nan(maxIter, 1);
end
history.restart = false(maxIter, 1);
end

function history = trimHistory(history, count)
fields = fieldnames(history);
for j = 1:numel(fields)
    history.(fields{j}) = history.(fields{j})(1:count, :);
end
end
