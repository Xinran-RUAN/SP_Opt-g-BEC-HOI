function [rho, info] = PositiveConservativeProx( ...
    z, tau, V, potreg, mass, h, options)
%POSITIVECONSERVATIVEPROX Potential prox with positivity and fixed mass.
%
% Solves
%   min 0.5*||rho-z||_h^2 + tau*h*sum(V.*p_sigma(rho))
%   s.t. rho>=0, h*sum(rho)=mass.
% The linear branch delegates directly to the Bo Lin semismooth positive
% conservative projection.

if nargin < 7 || isempty(options)
    options = struct();
end
options = applyDefaults(options);
wasRow = isrow(z);
z = z(:);
V = V(:);
potreg = normalizePotentialForProx(potreg);
validateInputs(z, tau, V, potreg, mass, h, options);

if strcmpi(potreg.prox_type, 'linear')
    linearSlope = potreg.dp_sigma(z);
    shifted = z - tau * V .* linearSlope;
    rho = src.constraints.ProjectPositiveConservative( ...
        shifted, mass, h, options.projection_tol);
    info = linearInfo(rho, z, tau * V .* linearSlope, mass, h);
    if wasRow
        rho = rho.';
    end
    return;
end

a = tau * V;
dpAtZero = potreg.dp_sigma(zeros(size(z)));
threshold = z - a .* dpAtZero;
lambdaHi = max(threshold);
span = max(1, max(threshold) - min(threshold));
lambdaLo = min(threshold) - span;
maxInnerIterations = 0;
[~, massLo, ~, innerIterations] = stateAtLambda( ...
    lambdaLo, z, a, potreg, h, options);
maxInnerIterations = max(maxInnerIterations, innerIterations);
bracketIterations = 0;
while massLo < mass && bracketIterations < options.bracket_max_iter
    span = 2 * span;
    lambdaLo = min(threshold) - span;
    [~, massLo, ~, innerIterations] = stateAtLambda( ...
        lambdaLo, z, a, potreg, h, options);
    maxInnerIterations = max(maxInnerIterations, innerIterations);
    bracketIterations = bracketIterations + 1;
end
if massLo < mass
    error('src:potential:PositiveConservativeProx:BracketFailure', ...
        'Unable to bracket the mass multiplier after %d expansions.', ...
        options.bracket_max_iter);
end

if isfinite(options.lambda_initial) ...
        && options.lambda_initial > lambdaLo ...
        && options.lambda_initial < lambdaHi
    lambda = options.lambda_initial;
else
    lambda = 0.5 * (lambdaLo + lambdaHi);
end
lambdaIterations = 0;
for iteration = 1:options.lambda_max_iter
    [~, currentMass, massDerivative, innerIterations] = stateAtLambda( ...
        lambda, z, a, potreg, h, options);
    maxInnerIterations = max(maxInnerIterations, innerIterations);
    lambdaIterations = iteration;
    massResidual = currentMass - mass;
    massError = abs(massResidual);
    if massError <= options.mass_tol
        break;
    end
    if massResidual > 0
        lambdaLo = lambda;
    else
        lambdaHi = lambda;
    end
    if isfinite(massDerivative) && massDerivative < 0
        lambdaNew = lambda - massResidual / massDerivative;
    else
        lambdaNew = NaN;
    end
    if ~isfinite(lambdaNew) || lambdaNew <= lambdaLo ...
            || lambdaNew >= lambdaHi
        lambdaNew = 0.5 * (lambdaLo + lambdaHi);
    end
    lambda = lambdaNew;
end

% Re-evaluate at the reported multiplier so rho, lambda, and KKT
% diagnostics refer to the same state even after the final Newton update.
[rho, currentMass, ~, innerIterations] = stateAtLambda( ...
    lambda, z, a, potreg, h, options);
maxInnerIterations = max(maxInnerIterations, innerIterations);
massError = abs(currentMass - mass);
info = genericInfo(rho, z, a, potreg, lambda, massError, ...
    lambdaIterations, maxInnerIterations, options);
if massError > 1e-12
    error('src:potential:PositiveConservativeProx:MassFailure', ...
        'Potential prox mass residual %.3e exceeds 1e-12.', massError);
end
if ~info.converged
    error('src:potential:PositiveConservativeProx:KKTFailure', ...
        ['Potential prox did not meet its mass/KKT tolerances: ' ...
        'mass %.3e, KKT %.3e.'], massError, info.kkt_residual);
end
if wasRow
    rho = rho.';
end
end

function [rho, currentMass, massDerivative, maxIterations] = ...
    stateAtLambda(lambda, z, a, potential, h, options)
dpAtZero = potential.dp_sigma(zeros(size(z)));
free = -z + a .* dpAtZero + lambda < 0;
rho = zeros(size(z));
maxIterations = 0;
if any(free)
    upper = z(free) - lambda;
    target = upper;
    lower = zeros(size(upper));
    afree = a(free);
    r = max(lower, upper - afree .* potential.dp_sigma(upper));
    converged = false(size(r));
    for iteration = 1:options.inner_max_iter
        residual = r - target + afree .* potential.dp_sigma(r);
        derivative = 1 + afree .* potential.d2p_sigma(r);
        newlyConverged = abs(residual) <= options.inner_tol;
        converged = converged | newlyConverged;
        if all(converged)
            maxIterations = iteration;
            break;
        end
        negative = residual < 0;
        lower(negative) = r(negative);
        upperBracket = ~negative;
        upper(upperBracket) = r(upperBracket);
        trial = r - residual ./ derivative;
        invalid = ~isfinite(trial) | trial <= lower | trial >= upper;
        trial(invalid) = 0.5 * (lower(invalid) + upper(invalid));
        r(~converged) = trial(~converged);
        maxIterations = iteration;
    end
    residual = r - (z(free) - lambda) ...
        + afree .* potential.dp_sigma(r);
    numericalFailureTolerance = max(10 * options.inner_tol, 1e-12);
    derivative = 1 + afree .* potential.d2p_sigma(r);
    numericalFailure = abs(residual) > numericalFailureTolerance;
    if any(numericalFailure)
        [maxResidual, worst] = max(abs(residual));
        worstRootCorrection = maxResidual / derivative(worst);
        error('src:potential:PositiveConservativeProx:InnerFailure', ...
            ['Safeguarded node Newton residual reached %.3e ' ...
            '(root correction %.3e, derivative %.3e, r %.3e, ' ...
            'a %.3e, target %.3e, bracket [%.3e, %.3e], iter %d).'], ...
            maxResidual, worstRootCorrection, derivative(worst), r(worst), ...
            afree(worst), target(worst), lower(worst), upper(worst), ...
            maxIterations);
    end
    rho(free) = r;
    derivative = 1 + afree .* potential.d2p_sigma(r);
    massDerivative = -h * sum(1 ./ derivative);
else
    massDerivative = 0;
end
currentMass = h * sum(rho);
end

function info = genericInfo(rho, z, a, potential, lambda, massError, ...
    lambdaIterations, maxInnerIterations, options)
free = -z + a .* potential.dp_sigma(zeros(size(z))) + lambda < 0;
active = ~free;
stationarity = rho - z + a .* potential.dp_sigma(rho) + lambda;
if any(free)
    freeResidual = abs(stationarity(free));
    freeStationarity = max(freeResidual);
    freeDerivative = 1 + a(free) .* potential.d2p_sigma(rho(free));
    freeRootErrors = freeResidual ./ freeDerivative;
    freeRootError = max(freeRootErrors);
else
    freeStationarity = 0;
    freeRootError = 0;
end
if any(active)
    dpAtZero = potential.dp_sigma(zeros(size(z)));
    activeDualViolation = max(max( ...
        z(active) - a(active) .* dpAtZero(active) - lambda, 0));
else
    activeDualViolation = 0;
end
info.lambda = lambda;
info.mass_error = massError;
info.active_count = nnz(active);
info.free_count = nnz(free);
info.lambda_iterations = lambdaIterations;
info.max_inner_iterations = maxInnerIterations;
info.free_stationarity = freeStationarity;
info.free_root_error = freeRootError;
info.active_dual_violation = activeDualViolation;
info.kkt_residual = max(freeStationarity, activeDualViolation);
kktTolerance = max(10 * options.inner_tol, 1e-12);
info.converged = massError <= options.mass_tol ...
    && info.kkt_residual <= kktTolerance;
info.backend = 'potential_positive_conservative';
end

function info = linearInfo(rho, z, linearShift, mass, h)
free = rho > 0;
active = ~free;
if any(free)
    lambda = mean(z(free) - linearShift(free) - rho(free));
else
    lambda = NaN;
end
stationarity = rho - z + linearShift + lambda;
if any(free)
    freeStationarity = max(abs(stationarity(free)));
else
    freeStationarity = 0;
end
if any(active)
    activeDualViolation = max(max(z(active) - linearShift(active) - lambda, 0));
else
    activeDualViolation = 0;
end
info.lambda = lambda;
info.mass_error = abs(h * sum(rho) - mass);
info.active_count = nnz(active);
info.free_count = nnz(free);
info.lambda_iterations = 0;
info.max_inner_iterations = 0;
info.free_stationarity = freeStationarity;
info.active_dual_violation = activeDualViolation;
info.kkt_residual = max(freeStationarity, activeDualViolation);
info.converged = info.mass_error <= 1e-12 ...
    && info.kkt_residual <= 1e-12;
info.backend = 'bolin_linear_potential';
end

function options = applyDefaults(options)
defaults.mass_tol = 1e-13;
defaults.lambda_max_iter = 50;
defaults.inner_tol = 1e-13;
defaults.inner_max_iter = 30;
defaults.bracket_max_iter = 60;
defaults.lambda_initial = NaN;
defaults.projection_tol = 1e-13;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
end

function validateInputs(z, tau, V, potreg, mass, h, options)
if isempty(z) || numel(V) ~= numel(z) || ~isreal(z) || ~isreal(V) ...
        || any(~isfinite(z)) || any(~isfinite(V))
    error('src:potential:PositiveConservativeProx:InvalidVectors', ...
        'z and V must be finite real vectors of equal nonzero length.');
end
if any(V < 0)
    error('src:potential:PositiveConservativeProx:NegativePotential', ...
        'Potential prox requires nonnegative V.');
end
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0 ...
        || ~isscalar(mass) || ~isfinite(mass) || mass <= 0 ...
        || ~isscalar(h) || ~isfinite(h) || h <= 0
    error('src:potential:PositiveConservativeProx:InvalidScalars', ...
        'tau, mass, and h must be positive finite scalars.');
end
required = {'p_sigma', 'dp_sigma', 'd2p_sigma', 'prox_type'};
if ~isstruct(potreg) || ~all(isfield(potreg, required)) ...
        || ~all(cellfun(@(field) isa(potreg.(field), 'function_handle'), ...
        required(1:3)))
    error('src:potential:PositiveConservativeProx:InvalidPotentialMap', ...
        'Potential prox requires p_sigma, dp_sigma, and d2p_sigma handles.');
end
positiveFields = {'mass_tol', 'lambda_max_iter', 'inner_tol', ...
    'inner_max_iter', 'bracket_max_iter', 'projection_tol'};
for j = 1:numel(positiveFields)
    value = options.(positiveFields{j});
    if ~isscalar(value) || ~isfinite(value) || value <= 0
        error('src:potential:PositiveConservativeProx:InvalidOptions', ...
            '%s must be a positive finite scalar.', positiveFields{j});
    end
end
end

function potential = normalizePotentialForProx(potential)
required = {'p_sigma', 'dp_sigma', 'd2p_sigma'};
if isstruct(potential) && all(isfield(potential, required))
    if ~isfield(potential, 'prox_type') || isempty(potential.prox_type)
        potential.prox_type = 'generic_convex';
    end
    return;
end
if ~isstruct(potential) || ~isfield(potential, 'name')
    error('src:potential:PositiveConservativeProx:InvalidPotentialMap', ...
        'Potential prox requires the generic handle interface.');
end
if strcmpi(potential.name, 'linear')
    potential = src.potential.MakeLinear();
    return;
end
if ~isfield(potential, 'sigma') || ~isscalar(potential.sigma) ...
        || ~isfinite(potential.sigma) || potential.sigma <= 0
    error('src:potential:PositiveConservativeProx:InvalidSigma', ...
        'Legacy sqrt prox input requires a positive finite sigma.');
end
sigma = potential.sigma;
potential.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
potential.dp_sigma = @(rho) rho ./ hypot(rho, sigma);
potential.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
potential.prox_type = 'generic_convex';
end
