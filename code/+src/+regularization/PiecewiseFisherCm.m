function regularization = PiecewiseFisherCm(epsilon, regularity)
%PIECEWISEFISHERCM Convexity-preserving C^m Fisher denominator, m=1,2.
%
% This construction acts on the density denominator r_epsilon(rho); it is
% not an amplitude regularization.  On 0<=rho<=epsilon,
%
% r = rho + epsilon - epsilon/(m+1)*(1-rho/epsilon)^(m+1),
%
% and r=rho+epsilon above the matching point.

if ~isscalar(epsilon) || ~isfinite(epsilon) || epsilon <= 0
    error('src:regularization:PiecewiseFisherCm:InvalidEpsilon', ...
        'epsilon must be a positive finite scalar.');
end
if ~isscalar(regularity) || ~ismember(regularity, [1, 2])
    error('src:regularization:PiecewiseFisherCm:InvalidRegularity', ...
        'Only regularity m=1 or m=2 is supported.');
end

regularization.r_epsilon = @(rho) evaluatePiecewise( ...
    rho, epsilon, regularity, 0);
regularization.dr_epsilon = @(rho) evaluatePiecewise( ...
    rho, epsilon, regularity, 1);
regularization.d2r_epsilon = @(rho) evaluatePiecewise( ...
    rho, epsilon, regularity, 2);
regularization.epsilon = epsilon;
regularization.name = sprintf('piecewise_c%d', regularity);
regularization.regularity = regularity;
regularization.label = sprintf([ ...
    'piecewise C^%d r_epsilon(rho), epsilon=%.6g'], ...
    regularity, epsilon);
regularization.source = 'PiecewiseFisherCm';
regularization.is_convex_preserving = true;
regularization = src.regularization.NormalizeFisher(regularization);
end

function value = evaluatePiecewise(rho, epsilon, regularity, order)
if ~isreal(rho) || any(~isfinite(rho(:)))
    error('src:regularization:PiecewiseFisherCm:InvalidDensity', ...
        'rho must be finite and real.');
end
inside = rho <= epsilon;
z = 1 - rho(inside) / epsilon;
switch order
    case 0
        value = rho + epsilon;
        value(inside) = rho(inside) + epsilon ...
            - epsilon / (regularity + 1) .* z .^ (regularity + 1);
    case 1
        value = ones(size(rho));
        value(inside) = 1 + z .^ regularity;
    case 2
        value = zeros(size(rho));
        value(inside) = -(regularity / epsilon) ...
            .* z .^ (regularity - 1);
    otherwise
        error('src:regularization:PiecewiseFisherCm:InternalOrder', ...
            'Unsupported derivative order.');
end
end
