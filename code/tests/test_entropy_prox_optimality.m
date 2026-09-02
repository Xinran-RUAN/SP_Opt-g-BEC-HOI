function stats = test_entropy_prox_optimality()
%TEST_ENTROPY_PROX_OPTIMALITY Common-lambda stationarity of the prox.

rng(37);
z = 0.05 * randn(48, 1);
h = 0.25;
mass = 2;
aValues = [1, 0.1, 0.01];
maximumResidual = 0;
for a = aValues
    [rho, info] = src.entropy.SimplexProx(z, a, mass, h, struct());
    assert(info.underflow_count == 0, ...
        'Optimality test unexpectedly entered an underflow case.');
    residual = rho - z + a * log(rho) + info.lambda;
    maximumResidual = max(maximumResidual, max(abs(residual)));
end
assert(maximumResidual <= 5e-13, ...
    'Entropy prox optimality residual is %.3e.', maximumResidual);
stats.maximum_optimality_error = maximumResidual;
fprintf('test_entropy_prox_optimality: max KKT %.3e\n', maximumResidual);
end
