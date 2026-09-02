function [rho, info] = SimplexProx(z, tauEta, mass, h, options)
%SIMPLEXPROX Prox of tauEta*H_h plus positivity and fixed mass.
%
% Solves rho-z+tauEta*log(rho)+lambda=0 without Lambert W. No density
% floor is imposed; exp(log(rho)) is allowed to underflow naturally.

if nargin < 5
    options = struct();
end
options = applyDefaults(options);
wasRow = isrow(z);
z = z(:);
if isempty(z) || any(~isfinite(z)) || ~isreal(z)
    error('src:entropy:SimplexProx:InvalidInput', ...
        'z must be a nonempty finite real vector.');
end
if ~isscalar(tauEta) || ~isfinite(tauEta) || tauEta <= 0
    error('src:entropy:SimplexProx:InvalidTauEta', ...
        'tau_eta must be a positive finite scalar.');
end
if ~isscalar(mass) || ~isfinite(mass) || mass <= 0 ...
        || ~isscalar(h) || ~isfinite(h) || h <= 0
    error('src:entropy:SimplexProx:InvalidMassOrSpacing', ...
        'mass and h must be positive finite scalars.');
end

a = tauEta;
meanDensity = mass / (h * numel(z));
phiMean = meanDensity + a * log(meanDensity);
lambdaLower = min(z) - phiMean;
lambdaUpper = max(z) - phiMean;
lambda = min(lambdaUpper, max(lambdaLower, mean(z) - phiMean));
innerMaximum = 0;

if lambdaLower == lambdaUpper
    [rho, logRho, innerIterations] = inversePhi(z - lambda, a, options);
    lambdaIterations = 0;
    innerMaximum = innerIterations;
else
    lambdaIterations = 0;
    for iteration = 1:options.lambda_max_iter
        lambdaIterations = iteration;
        [rho, ~, innerIterations] = inversePhi(z - lambda, a, options);
        innerMaximum = max(innerMaximum, innerIterations);
        massValue = h * sum(rho);
        massDifference = massValue - mass;
        if abs(massDifference) <= options.mass_tol
            break;
        end

        if massDifference > 0
            lambdaLower = lambda;
        else
            lambdaUpper = lambda;
        end
        derivative = -h * sum(rho ./ (rho + a));
        if isfinite(derivative) && derivative < 0
            newtonLambda = lambda - massDifference / derivative;
        else
            newtonLambda = NaN;
        end
        if ~isfinite(newtonLambda) || newtonLambda <= lambdaLower ...
                || newtonLambda >= lambdaUpper
            lambda = lambdaLower + 0.5 * (lambdaUpper - lambdaLower);
        else
            lambda = newtonLambda;
        end
    end
    [rho, logRho, innerIterations] = inversePhi(z - lambda, a, options);
    innerMaximum = max(innerMaximum, innerIterations);
end

info.mass_error = abs(h * sum(rho) - mass);
info.lambda = lambda;
info.lambda_iterations = lambdaIterations;
info.inner_max_iterations = innerMaximum;
info.min_rho = min(rho);
info.min_log_rho = min(logRho);
info.underflow_count = nnz(rho == 0 & isfinite(logRho));
info.converged = info.mass_error <= options.mass_tol;
if wasRow
    rho = rho.';
end
end

function [rho, u, maximumIterations] = inversePhi(b, a, options)
% Brackets follow from F(0)=1-b and F((b-1)/a)<=0 when b<=1.
lower = zeros(size(b));
upper = zeros(size(b));
large = b > 1;
lower(large) = 0;
upper(large) = log(b(large));
lower(~large) = max(-realmax, (b(~large) - 1) / a);
upper(~large) = 0;

u = zeros(size(b));
recommendedLog = b > a;
u(recommendedLog) = log(max(b(recommendedLog), realmin));
u(~recommendedLog) = max(-realmax, min(realmax, b(~recommendedLog) / a));
u = min(upper, max(lower, u));
maximumIterations = 0;

for iteration = 1:options.inner_newton_max_iter
    exponential = exp(u);
    residual = exponential + a * u - b;
    scale = max(1, abs(b) + exponential + abs(a * u));
    converged = abs(residual) <= options.inner_tol .* scale;
    if all(converged)
        maximumIterations = iteration - 1;
        break;
    end

    positive = residual > 0;
    upper(positive) = u(positive);
    lower(~positive) = u(~positive);
    newton = u - residual ./ (exponential + a);
    unsafe = ~isfinite(newton) | newton <= lower | newton >= upper;
    newton(unsafe) = lower(unsafe) ...
        + 0.5 * (upper(unsafe) - lower(unsafe));
    stepConverged = abs(newton - u) <= options.inner_tol ...
        .* max(1, abs(u));
    u(~converged) = newton(~converged);
    maximumIterations = iteration;
    if all(converged | stepConverged)
        break;
    end
end
rho = exp(u);
end

function options = applyDefaults(options)
defaults.mass_tol = 1e-14;
defaults.lambda_max_iter = 50;
defaults.inner_newton_max_iter = 30;
defaults.inner_tol = 1e-14;
names = fieldnames(defaults);
for j = 1:numel(names)
    if ~isfield(options, names{j}) || isempty(options.(names{j}))
        options.(names{j}) = defaults.(names{j});
    end
end
end
