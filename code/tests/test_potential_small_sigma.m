function stats = test_potential_small_sigma()
%TEST_POTENTIAL_SMALL_SIGMA Stable value/derivatives down to sigma=1e-12.

rho = [0; 1e-16; 1e-14; 1e-12; 1e-10; 1e-8; 1e-6; 1e-3];
sigmaList = [1e-3, 1e-6, 1e-9, 1e-12];
maxDirectRelativeError = 0;
maxDerivativeBoundViolation = 0;
for j = 1:numel(sigmaList)
    sigma = sigmaList(j);
    expression = src.potential.Expression('sqrt_power');
    p = expression.value(rho, sigma);
    dp = expression.first_derivative(rho, sigma);
    d2p = expression.second_derivative(rho, sigma);

    assert(all(isfinite(p)) && all(isfinite(dp)) && all(isfinite(d2p)), ...
        'Potential evaluation is nonfinite for sigma %.1e.', sigma);
    assert(all(p >= 0) && all(p <= rho), ...
        'Potential value bounds failed for sigma %.1e.', sigma);
    boundViolation = max([max(-dp), max(dp - 1), max(-d2p)]);
    maxDerivativeBoundViolation = max( ...
        maxDerivativeBoundViolation, boundViolation);
    assert(boundViolation <= 10 * eps, ...
        'Potential derivative bounds failed for sigma %.1e.', sigma);
    assert(dp(1) == 0 && abs(d2p(1) - 1 / sigma) <= 10 * eps / sigma, ...
        'Origin derivatives failed for sigma %.1e.', sigma);

    % Compare the subtractive formula only where cancellation is benign.
    direct = hypot(rho, sigma) - sigma;
    safe = rho >= 1e-2 * sigma & rho > 0;
    if any(safe)
        relativeError = max(abs(p(safe) - direct(safe)) ...
            ./ max(p(safe), realmin));
        maxDirectRelativeError = max(maxDirectRelativeError, relativeError);
        assert(relativeError <= 1e-10, ...
            'Stable/direct value disagreement %.3e for sigma %.1e.', ...
            relativeError, sigma);
    end
end

epsilon = 1e-3;
for power = 1:4
    [potreg, info] = src.potential.Validate( ...
        struct('name', 'sqrt_power', 'power', power), epsilon);
    assert(potreg.sigma == epsilon ^ power && info.power == power, ...
        'sqrt_power scale mapping failed for p=%d.', power);
end

stats.sigma = sigmaList;
stats.max_direct_relative_error = maxDirectRelativeError;
stats.max_derivative_bound_violation = maxDerivativeBoundViolation;
stats.all_finite = true;
fprintf(['test_potential_small_sigma: finite PASS, direct rel %.3e, ' ...
    'bound violation %.3e\n'], maxDirectRelativeError, ...
    maxDerivativeBoundViolation);
end
