function result = PolishEntropyNewton(rho0, problem, options)
%POLISHENTROPYNEWTON Matrix-free equality-constrained interior Newton polish.

options = applyDefaults(options);
startTime = tic;
rho = rho0(:);
energy = src.discretization.ps.Energy(rho, problem);
gradient = src.discretization.ps.Gradient(rho, problem);
state = src.solvers.EvaluateState(rho, problem, options, energy, gradient);
initialState = state;
history = allocateHistory(options.max_iter);
iterationCount = 0;
attempted = false;
failed = false;
failureMessage = '';
status = 'not_attempted';

if any(rho <= 0)
    failed = true;
    status = 'entropy_newton_underflow';
    failureMessage = ['Entropy Newton skipped because at least one density ' ...
        'entry underflowed to zero.'];
    lambda = NaN;
    kkt = src.entropy.KKTResidual(rho, gradient, problem);
    converged = false;
    result = packageResult();
    return;
end
if state.pg_residual > options.entry_pg_tol
    status = 'entry_rejected';
    failureMessage = sprintf( ...
        'Entropy Newton entry rejected: PG residual %.3e exceeds %.3e.', ...
        state.pg_residual, options.entry_pg_tol);
    lambda = -mean(gradient + problem.entropy.eta * log(rho));
    kkt = src.entropy.KKTResidual(rho, gradient, problem, lambda);
    converged = false;
    result = packageResult();
    return;
end

attempted = true;
lambda = -mean(gradient + problem.entropy.eta * log(rho));
kkt = src.entropy.KKTResidual(rho, gradient, problem, lambda);
converged = kkt.kkt_residual <= options.pg_tol;
status = 'running';

for iteration = 1:options.max_iter
    if converged
        status = 'converged';
        break;
    end
    stationarity = gradient + problem.entropy.eta * log(rho) + lambda;
    massResidual = src.constraints.Mass(rho, problem.grid.h) - problem.mass;
    rightHandSide = -[stationarity; massResidual];
    operator = @(z) kktAction(z, rho, problem);
    numberOfUnknowns = problem.grid.N + 1;
    restart = min(50, numberOfUnknowns);
    maxOuter = min(numberOfUnknowns, ...
        max(1, ceil(options.gmres_maxit / restart)));
    [correction, gmresFlag, gmresRelativeResidual, gmresIteration] = gmres( ...
        operator, rightHandSide, restart, options.gmres_tol, maxOuter);
    if gmresIteration(1) == 0
        gmresIterations = gmresIteration(2);
    else
        gmresIterations = (gmresIteration(1) - 1) * restart ...
            + gmresIteration(2);
    end
    if gmresFlag ~= 0 && gmresRelativeResidual > options.gmres_failure_tol
        failed = true;
        status = 'gmres_failed';
        failureMessage = sprintf( ...
            'Entropy Newton GMRES failed: flag=%d, relative residual=%.3e.', ...
            gmresFlag, gmresRelativeResidual);
        break;
    end

    direction = correction(1:end-1);
    lambdaDirection = correction(end);
    % Enforce the scalar Newton mass equation to working precision even
    % when the matrix-free GMRES solve is deliberately inexact.
    massDirectionResidual = problem.grid.h * sum(direction) + massResidual;
    direction = direction - massDirectionResidual ...
        / (problem.grid.h * problem.grid.N);
    negative = direction < 0;
    if any(negative)
        maximumPositiveStep = min(-rho(negative) ./ direction(negative));
    else
        maximumPositiveStep = inf;
    end
    step = min(1, 0.99 * maximumPositiveStep);
    accepted = false;
    backtracks = 0;
    currentResidual = kkt.kkt_residual;

    while step >= options.min_step && backtracks <= options.max_backtracks
        trialRho = rho + step * direction;
        trialLambda = lambda + step * lambdaDirection;
        if all(trialRho > 0)
            trialMassError = abs(src.constraints.Mass( ...
                trialRho, problem.grid.h) - problem.mass);
            if trialMassError <= options.mass_tol
                trialEnergy = src.discretization.ps.Energy(trialRho, problem);
                trialGradient = src.discretization.ps.Gradient(trialRho, problem);
                trialKkt = src.entropy.KKTResidual( ...
                    trialRho, trialGradient, problem, trialLambda);
                residualBound = (1 - options.residual_c1 * step) ...
                    * currentResidual;
                if trialKkt.kkt_residual <= residualBound ...
                        || trialKkt.kkt_residual <= options.pg_tol
                    accepted = true;
                    break;
                end
            end
        end
        step = options.backtrack * step;
        backtracks = backtracks + 1;
    end

    if ~accepted
        failed = true;
        status = 'globalization_failed';
        failureMessage = 'Entropy Newton residual globalization failed.';
        break;
    end

    rho = trialRho;
    lambda = trialLambda;
    energy = trialEnergy;
    gradient = trialGradient;
    kkt = trialKkt;
    state = src.solvers.EvaluateState(rho, problem, options, energy, gradient);
    iterationCount = iteration;
    history.iteration(iteration) = iteration;
    history.elapsed_time(iteration) = toc(startTime);
    history.physical_energy(iteration) = state.physical_energy;
    history.entropy_value(iteration) = state.entropy_value;
    history.augmented_energy(iteration) = state.augmented_energy;
    history.pg_residual(iteration) = state.pg_residual;
    history.kkt_residual(iteration) = kkt.kkt_residual;
    history.mass_error(iteration) = state.mass_error;
    history.min_density(iteration) = state.min_density;
    history.gmres_iterations(iteration) = gmresIterations;
    history.gmres_flag(iteration) = gmresFlag;
    history.gmres_relative_residual(iteration) = gmresRelativeResidual;
    history.step(iteration) = step;
    history.backtracks(iteration) = backtracks;

    if options.display
        fprintf('EN %3d  %.3e  %.3e  %4d  %.3e\n', ...
            iteration, state.pg_residual, kkt.kkt_residual, ...
            gmresIterations, step);
    end
    converged = kkt.kkt_residual <= options.pg_tol;
    if converged
        status = 'converged';
    end
end
if ~failed && ~converged && strcmp(status, 'running')
    status = 'max_iterations';
end
result = packageResult();

    function output = packageResult()
        output.rho = rho;
        output.energy = energy;
        output.history = trimHistory(history, iterationCount);
        output.polish_attempted = attempted;
        output.polish_converged = converged;
        output.failed = failed;
        output.failure_message = failureMessage;
        output.status = status;
        output.iterations = iterationCount;
        output.elapsed_time = toc(startTime);
        output.pg_residual_before = initialState.pg_residual;
        output.pg_residual = state.pg_residual;
        output.kkt_residual = kkt.kkt_residual;
        output.mass_error = state.mass_error;
        output.min_density = state.min_density;
        output.lambda = lambda;
        output.energy_change = energy - initialState.energy;
        output.state_change = sqrt(problem.grid.h ...
            * sum((rho - rho0(:)) .^ 2));
    end
end

function value = kktAction(z, rho, problem)
direction = z(1:end-1);
hessianDirection = src.solvers.HessianAction( ...
    rho, direction, problem, problem.plan, problem.regularization);
hessianDirection = hessianDirection ...
    + problem.entropy.eta * (direction ./ rho);
value = [hessianDirection + z(end); ...
    problem.grid.h * sum(direction)];
end

function options = applyDefaults(options)
if nargin < 1 || isempty(options)
    options = struct();
end
defaults.entry_pg_tol = 1e-5;
defaults.pg_tol = 1e-12;
defaults.max_iter = 50;
defaults.gmres_tol = 1e-10;
defaults.gmres_maxit = 200;
defaults.gmres_failure_tol = 0.5;
defaults.residual_step = 1;
defaults.projection_name = 'simplex';
defaults.projection_tol = 1e-13;
defaults.residual_c1 = 1e-4;
defaults.backtrack = 0.5;
defaults.min_step = 1e-12;
defaults.max_backtracks = 30;
defaults.mass_tol = 1e-10;
defaults.display = true;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
end

function history = allocateHistory(maxIter)
fields = {'iteration', 'elapsed_time', 'physical_energy', ...
    'entropy_value', 'augmented_energy', 'pg_residual', 'kkt_residual', ...
    'mass_error', 'min_density', 'gmres_iterations', 'gmres_flag', ...
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
