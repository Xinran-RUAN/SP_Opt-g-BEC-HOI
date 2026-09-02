function [r, dr, d2r, info] = EvaluateDenominator(rho, reg)
%EVALUATEDENOMINATOR Legacy named adapter for r_epsilon derivatives.
%
% The active definitions have physical domain rho >= 0. Values that are
% negative only at roundoff scale are clipped to zero; material negative
% values are rejected rather than extending the regularization.

[reg, info] = src.regularization.Validate(reg);
if ~isreal(rho) || any(~isfinite(rho(:)))
    error('src:regularization:EvaluateDenominator:InvalidDensity', ...
        'rho must be finite and real.');
end
roundoff = 100 * eps(max(1, max(abs(rho(:)))));
if any(rho(:) < -roundoff)
    error('src:regularization:EvaluateDenominator:NegativeDensity', ...
        'Active regularizations are defined only for rho >= 0.');
end
r = max(rho, 0);

switch reg.name
    case 'shift_smooth'
        denominator = r + reg.epsilon;
        derivative = ones(size(r));
        secondDerivative = zeros(size(r));
    case {'piecewise_c1', 'piecewise_c2'}
        m = sscanf(reg.name, 'piecewise_c%d');
        fisher = src.regularization.PiecewiseFisherCm(reg.epsilon, m);
        [denominator, derivative, secondDerivative] = ...
            src.regularization.EvaluateFisher(r, fisher);
    otherwise
        m = sscanf(reg.name, 'piecewise_c%d');
        d = reg.transition_width;
        inside = (r <= d);
        z = 1 - r(inside) / d;

        denominator = reg.epsilon + r + d / (m + 1);
        derivative = ones(size(r));
        secondDerivative = zeros(size(r));
        denominator(inside) = reg.epsilon + r(inside) + ...
            d / (m + 1) .* (1 - z .^ (m + 1));
        derivative(inside) = 1 + z .^ m;
        secondDerivative(inside) = -(m / d) .* z .^ (m - 1);
end

if any(denominator(:) <= 0) || any(secondDerivative(:) ...
        > 100 * eps(max(1, max(abs(secondDerivative(:))))))
    error('src:regularization:EvaluateDenominator:InvariantFailure', ...
        'The denominator must be positive and concave.');
end
r = denominator;
dr = derivative;
d2r = secondDerivative;
end
