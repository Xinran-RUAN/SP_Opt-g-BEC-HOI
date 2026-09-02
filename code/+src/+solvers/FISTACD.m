function [rho, diagnostics, history] = FISTACD(problem, rho0, solver)
%FISTACD Optional accelerated comparison method with robust restarts.
%
% Feasibility and monotone restarts alter the standard accelerated method;
% this implementation does not claim all unrestarted FISTA rate theory.

startTime = tic;
rho = rho0(:);
rhoPrevious = rho;
energy = src.discretization.ps.Energy(rho, problem);
gradient = src.discretization.ps.Gradient(rho, problem);
[smoothEnergy, smoothGradient] = splittingQuantities( ...
    rho, energy, gradient, problem, solver);
state = src.solvers.EvaluateState(rho, problem, solver, energy, gradient);
objective = state.augmented_energy;
L = solver.L0;
history = allocateHistory(solver.max_iter);
switchEnabled = solver.switch.enabled;
stopAtEnergyHandoff = switchEnabled ...
    && solver.switch.stop_at_energy_handoff;
if stopAtEnergyHandoff
    mainIterationLimit = min(solver.max_iter, solver.switch.max_main_iter);
else
    mainIterationLimit = solver.max_iter;
end
fullPgResidual = fullResidual(rho, gradient, state, problem, solver);
if switchEnabled
    converged = fullPgResidual <= solver.final_pg_tol;
else
    converged = state.pg_residual <= solver.pg_tol;
end
restartCount = 0;
relativeEnergyChange = NaN;
iterationCount = 0;
windowEnergySpan = NaN;
plateauCounter = 0;
energyPlateauDetected = false;
requestPolish = false;
handoffCaptured = false;
handoffIteration = NaN;
handoffTime = NaN;
handoffEnergy = NaN;
handoffPgResidual = NaN;
handoffRho = [];
if converged && switchEnabled
    stopReason = 'final_pg_tolerance';
else
    stopReason = 'maximum_iterations';
end

for iteration = 1:mainIterationLimit
    if converged
        break;
    end
    alpha = (iteration - 1) / (iteration + solver.a);
    extrapolated = rho + alpha * (rho - rhoPrevious);
    restarted = false;

    if usesPotentialProx(problem, solver) && any(extrapolated < 0)
        extrapolated = rho;
        restarted = true;
    elseif any(extrapolated < -solver.feasibility_tol)
        extrapolated = rho;
        restarted = true;
    elseif any(extrapolated < 0)
        extrapolated = projectFeasible(max(extrapolated, 0), problem, solver);
    end

    [smoothEnergyY, smoothGradientY] = splittingQuantities( ...
        extrapolated, [], [], problem, solver);
    [candidate, candidateSmoothEnergy, L, backtracks, proxInfo] = ...
        backtrackingStep(extrapolated, smoothEnergyY, smoothGradientY, ...
        L, problem, solver);
    candidateEnergy = src.discretization.ps.Energy(candidate, problem);

    candidateEntropy = entropyValue(candidate, problem);
    candidateObjective = candidateEnergy + entropyEta(problem) * candidateEntropy;
    monotoneTolerance = 100 * eps(max([1, abs(objective), abs(candidateObjective)]));
    if candidateObjective > objective + monotoneTolerance
        extrapolated = rho;
        smoothEnergyY = smoothEnergy;
        smoothGradientY = smoothGradient;
        [candidate, candidateSmoothEnergy, L, restartBacktracks, proxInfo] = ...
            backtrackingStep(extrapolated, smoothEnergyY, smoothGradientY, ...
            L, problem, solver);
        backtracks = backtracks + restartBacktracks;
        candidateEnergy = src.discretization.ps.Energy(candidate, problem);
        candidateEntropy = entropyValue(candidate, problem);
        candidateObjective = candidateEnergy ...
            + entropyEta(problem) * candidateEntropy;
        restarted = true;
    end

    candidateGradient = src.discretization.ps.Gradient(candidate, problem);
    [candidateSmoothEnergy, candidateSmoothGradient] = ...
        splittingQuantities(candidate, candidateEnergy, candidateGradient, ...
        problem, solver, candidateSmoothEnergy);
    residualChecked = iteration == 1 ...
        || mod(iteration, solver.residual_check_interval) == 0;
    candidateState = src.solvers.EvaluateState(candidate, problem, solver, ...
        candidateEnergy, candidateGradient, residualChecked);
    if ~residualChecked
        candidateState.pg_residual = state.pg_residual;
    end
    if switchEnabled
        candidateFullPgResidual = fullResidual( ...
            candidate, candidateGradient, candidateState, problem, solver);
    else
        candidateFullPgResidual = NaN;
    end
    relativeEnergyChange = abs(candidateObjective - objective) ...
        / max(1, abs(objective));
    if restarted
        restartCount = restartCount + 1;
    end

    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.energy(iteration) = candidateEnergy;
    history.physical_energy(iteration) = candidateEnergy;
    history.entropy_value(iteration) = candidateEntropy;
    history.augmented_energy(iteration) = candidateObjective;
    history.pg_residual(iteration) = candidateState.pg_residual;
    history.full_pg_residual(iteration) = candidateFullPgResidual;
    history.relative_energy_change(iteration) = relativeEnergyChange;
    history.mass_error(iteration) = candidateState.mass_error;
    history.min_density(iteration) = candidateState.min_density;
    history.accepted_L(iteration) = L;
    history.L(iteration) = L;
    history.tau(iteration) = 1 / L;
    history.backtracks(iteration) = backtracks;
    history.prox_lambda_iterations(iteration) = ...
        proxIteration(proxInfo, 'lambda');
    history.prox_inner_iterations(iteration) = ...
        proxIteration(proxInfo, 'inner');
    history.residual_checked(iteration) = residualChecked;
    history.restart(iteration) = restarted;

    if switchEnabled && iteration >= solver.switch.energy_window
        firstWindowIndex = iteration - solver.switch.energy_window + 1;
        energyWindow = history.augmented_energy(firstWindowIndex:iteration);
        windowEnergySpan = (max(energyWindow) - min(energyWindow)) ...
            / max(1, abs(candidateObjective));
    else
        windowEnergySpan = NaN;
    end
    energyPlateau = switchEnabled ...
        && iteration >= solver.switch.min_iter ...
        && isfinite(windowEnergySpan) ...
        && windowEnergySpan <= solver.switch.energy_tol;
    energyPlateauDetected = energyPlateauDetected || energyPlateau;
    localEnough = switchEnabled ...
        && candidateFullPgResidual <= solver.switch.pg_entry_tol;
    if energyPlateau && localEnough
        plateauCounter = plateauCounter + 1;
    else
        plateauCounter = 0;
    end
    history.window_energy_span(iteration) = windowEnergySpan;
    history.energy_plateau_detected(iteration) = energyPlateau;
    history.plateau_counter(iteration) = plateauCounter;
    handoffDetected = plateauCounter >= solver.switch.consecutive_windows;
    history.handoff_detected(iteration) = handoffDetected;
    if handoffDetected && ~handoffCaptured
        handoffCaptured = true;
        handoffIteration = iteration;
        handoffTime = history.elapsed_time(iteration);
        handoffEnergy = candidateObjective;
        handoffPgResidual = candidateFullPgResidual;
        if solver.history.capture_handoff_state
            handoffRho = candidate;
        end
    end

    oldRho = rho;
    rho = candidate;
    energy = candidateEnergy;
    objective = candidateObjective;
    gradient = candidateGradient;
    smoothEnergy = candidateSmoothEnergy;
    smoothGradient = candidateSmoothGradient;
    state = candidateState;
    fullPgResidual = candidateFullPgResidual;
    if restarted
        rhoPrevious = rho;
    else
        rhoPrevious = oldRho;
    end

    displayInterval = solver.display_interval;
    if switchEnabled
        displayInterval = solver.display_every;
    end
    if solver.display && (iteration == 1 || mod(iteration, displayInterval) == 0)
        if switchEnabled
            fprintf('%7d  %.15e  %.3e  %.3e  %.3e  %d\n', ...
                iteration, energy, fullPgResidual, windowEnergySpan, L, restarted);
        else
        fprintf('%7d  %.15e  %.3e  %.3e  %.3e  %d\n', ...
            iteration, energy, state.pg_residual, relativeEnergyChange, L, restarted);
        end
    end
    if switchEnabled
        if fullPgResidual <= solver.final_pg_tol
            converged = true;
            stopReason = 'final_pg_tolerance';
            break;
        elseif stopAtEnergyHandoff && handoffDetected
            requestPolish = true;
            stopReason = 'switch_to_kkt_polish';
            break;
        end
    elseif residualChecked && state.pg_residual <= solver.pg_tol
        converged = true;
        stopReason = 'pg_tolerance';
    end
end

if stopAtEnergyHandoff && ~converged && ~requestPolish ...
        && iterationCount >= mainIterationLimit
    if fullPgResidual <= solver.switch.forced_pg_tol
        requestPolish = true;
        stopReason = 'forced_local_polish';
    else
        stopReason = 'main_solver_not_in_local_regime';
    end
end

history = trimHistory(history, iterationCount);
if iterationCount > 0 && ~history.residual_checked(end)
    state = src.solvers.EvaluateState( ...
        rho, problem, solver, energy, gradient, true);
    history.pg_residual(end) = state.pg_residual;
    history.residual_checked(end) = true;
end
diagnostics.iterations = iterationCount;
diagnostics.energy = energy;
diagnostics.physical_energy = state.physical_energy;
diagnostics.entropy_value = state.entropy_value;
diagnostics.augmented_energy = state.augmented_energy;
diagnostics.pg_residual = state.pg_residual;
diagnostics.main_stop_reason = stopReason;
diagnostics.stop_reason = stopReason;
diagnostics.request_polish = requestPolish;
diagnostics.energy_window_span = windowEnergySpan;
diagnostics.main_energy_window_span = windowEnergySpan;
diagnostics.energy_plateau_detected = energyPlateauDetected;
diagnostics.plateau_counter = plateauCounter;
diagnostics.kkt_residual = state.kkt_residual;
diagnostics.relative_energy_change = relativeEnergyChange;
diagnostics.mass_error = state.mass_error;
diagnostics.min_density = state.min_density;
diagnostics.accepted_L = L;
diagnostics.accepted_L_final = L;
diagnostics.final_tau = 1 / L;
if isempty(history.accepted_L)
    diagnostics.accepted_L_max = L;
    diagnostics.mean_backtracks = 0;
    diagnostics.max_backtracks = 0;
else
    diagnostics.accepted_L_max = max(history.accepted_L);
    diagnostics.mean_backtracks = mean(history.backtracks);
    diagnostics.max_backtracks = max(history.backtracks);
end
if isempty(history.prox_lambda_iterations)
    diagnostics.mean_prox_lambda_iterations = 0;
    diagnostics.max_prox_lambda_iterations = 0;
    diagnostics.mean_prox_inner_iterations = 0;
    diagnostics.max_prox_inner_iterations = 0;
else
    diagnostics.mean_prox_lambda_iterations = ...
        mean(history.prox_lambda_iterations);
    diagnostics.max_prox_lambda_iterations = ...
        max(history.prox_lambda_iterations);
    diagnostics.mean_prox_inner_iterations = ...
        mean(history.prox_inner_iterations);
    diagnostics.max_prox_inner_iterations = ...
        max(history.prox_inner_iterations);
end
if entropyActive(problem)
    diagnostics.composite_pg_residual = state.pg_residual;
    diagnostics.full_pg_residual = NaN;
else
    diagnostics.composite_pg_residual = ...
        src.solvers.CompositeGradientMapping(rho, problem, solver);
    diagnostics.full_pg_residual = src.solvers.FullGradientMapping( ...
        rho, gradient, problem, solver);
end
if switchEnabled
    diagnostics.full_pg_residual = fullPgResidual;
end
diagnostics.converged = converged;
diagnostics.failed = false;
diagnostics.failure_message = '';
diagnostics.restarts = restartCount;
diagnostics.restart_count = restartCount;
diagnostics.elapsed_time = toc(startTime);
diagnostics.energy_plateau = energyPlateauDetected;
diagnostics.stop_at_energy_handoff = stopAtEnergyHandoff;
diagnostics.handoff_detected = handoffCaptured;
diagnostics.handoff_iteration = handoffIteration;
diagnostics.handoff_time = handoffTime;
diagnostics.handoff_energy = handoffEnergy;
diagnostics.handoff_pg_residual = handoffPgResidual;
diagnostics.handoff_rho = handoffRho;
end

function rho = projectFeasible(z, problem, solver)
switch lower(char(solver.projection_name))
    case 'simplex'
        rho = src.constraints.ProjectSimplex(z, problem.mass, problem.grid.h);
    case 'semismooth'
        rho = src.constraints.ProjectPositiveConservative( ...
            z, problem.mass, problem.grid.h, solver.projection_tol);
end
end

function [candidate, candidateSmoothEnergy, acceptedL, backtracks, proxInfo] = ...
    backtrackingStep(y, smoothEnergyY, smoothGradientY, initialL, ...
    problem, solver)
acceptedL = initialL;
for trial = 1:solver.max_backtracks
    [candidate, proxInfo] = src.solvers.ProxStep( ...
        y, smoothGradientY, acceptedL, problem, solver);
    candidateSmoothEnergy = splittingEnergy(candidate, problem, solver);
    difference = candidate - y;
    majorant = smoothEnergyY ...
        + problem.grid.h * sum(smoothGradientY .* difference) ...
        + 0.5 * acceptedL * problem.grid.h * sum(difference .^ 2);
    roundoff = 100 * eps(max( ...
        [1, abs(smoothEnergyY), abs(candidateSmoothEnergy)]));
    if candidateSmoothEnergy <= majorant + roundoff
        backtracks = trial - 1;
        return;
    end
    acceptedL = solver.backtrack_factor * acceptedL;
    if ~isfinite(acceptedL)
        break;
    end
end
error('src:solvers:FISTACD:BacktrackingFailure', ...
    'Backtracking failed to establish the h-inner-product majorization.');
end

function history = allocateHistory(maxIter)
fields = {'iteration', 'elapsed_time', 'energy', 'physical_energy', ...
    'entropy_value', 'augmented_energy', 'pg_residual', ...
    'full_pg_residual', 'window_energy_span', 'plateau_counter', ...
    'relative_energy_change', 'mass_error', 'min_density', 'accepted_L', ...
    'L', 'tau', 'backtracks', 'prox_lambda_iterations', ...
    'prox_inner_iterations'};
for j = 1:numel(fields)
    history.(fields{j}) = nan(maxIter, 1);
end
history.restart = false(maxIter, 1);
history.residual_checked = false(maxIter, 1);
history.energy_plateau_detected = false(maxIter, 1);
history.handoff_detected = false(maxIter, 1);
end

function residual = fullResidual(rho, gradient, state, problem, solver)
if entropyActive(problem)
    residual = state.pg_residual;
else
    residual = src.solvers.FullGradientMapping( ...
        rho, gradient, problem, solver);
end
end

function history = trimHistory(history, count)
fields = fieldnames(history);
for j = 1:numel(fields)
    history.(fields{j}) = history.(fields{j})(1:count, :);
end
end

function eta = entropyEta(problem)
eta = 0;
if isfield(problem, 'entropy') && isfield(problem.entropy, 'enabled') ...
        && problem.entropy.enabled && isfield(problem.entropy, 'eta') ...
        && problem.entropy.eta > 0
    eta = problem.entropy.eta;
end
end

function value = entropyValue(rho, problem)
value = src.entropy.Value(rho, problem.grid);
end

function [energy, gradient] = splittingQuantities( ...
    rho, fullEnergy, fullGradient, problem, solver, knownSmoothEnergy)
if nargin < 6
    knownSmoothEnergy = [];
end
if usesPotentialProx(problem, solver)
    if isempty(knownSmoothEnergy)
        energy = src.discretization.ps.SmoothEnergy(rho, problem);
    else
        energy = knownSmoothEnergy;
    end
    gradient = src.discretization.ps.SmoothGradient(rho, problem);
else
    if isempty(fullEnergy)
        energy = src.discretization.ps.Energy(rho, problem);
    else
        energy = fullEnergy;
    end
    if isempty(fullGradient)
        gradient = src.discretization.ps.Gradient(rho, problem);
    else
        gradient = fullGradient;
    end
end
end

function energy = splittingEnergy(rho, problem, solver)
if usesPotentialProx(problem, solver)
    energy = src.discretization.ps.SmoothEnergy(rho, problem);
else
    energy = src.discretization.ps.Energy(rho, problem);
end
end

function active = usesPotentialProx(problem, solver)
active = ~entropyActive(problem) && isfield(solver, 'splitting') ...
    && strcmpi(solver.splitting, 'potential_prox');
end

function active = entropyActive(problem)
active = isfield(problem, 'entropy') ...
    && isfield(problem.entropy, 'enabled') && problem.entropy.enabled ...
    && isfield(problem.entropy, 'eta') && problem.entropy.eta > 0;
end

function count = proxIteration(info, kind)
switch kind
    case 'lambda'
        field = 'lambda_iterations';
    case 'inner'
        if isfield(info, 'max_inner_iterations')
            field = 'max_inner_iterations';
        else
            field = 'inner_max_iterations';
        end
end
if isfield(info, field)
    count = info.(field);
else
    count = 0;
end
end
