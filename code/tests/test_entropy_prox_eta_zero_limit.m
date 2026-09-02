function stats = test_entropy_prox_eta_zero_limit()
%TEST_ENTROPY_PROX_ETA_ZERO_LIMIT Convergence to the ordinary simplex prox.

rng(41);
z = 0.1 * randn(64, 1);
h = 0.2;
mass = 1;
etaValues = 10 .^ (-(2:10));
simplex = src.constraints.ProjectSimplex(z, mass, h);
errors = zeros(size(etaValues));
for j = 1:numel(etaValues)
    entropyState = src.entropy.SimplexProx( ...
        z, etaValues(j), mass, h, struct());
    errors(j) = sqrt(h * sum((entropyState - simplex) .^ 2));
end
assert(errors(end) < 1e-3 * errors(1), ...
    'Entropy prox did not approach the simplex projection.');
assert(min(diff(errors)) <= 0 && errors(end) <= min(errors(1:3)), ...
    'Entropy prox errors do not show a decreasing trend.');

problem.mass = mass;
problem.grid.h = h;
problem.entropy.enabled = true;
problem.entropy.eta = 0;
solver.projection_name = 'simplex';
solver.projection_tol = 1e-13;
etaZeroState = src.solvers.ApplyProx(z, 1, problem, solver);
stats.eta_zero_bitwise_equal = isequal(etaZeroState, simplex);
assert(stats.eta_zero_bitwise_equal, ...
    'eta=0 did not use the exact existing simplex path.');
stats.eta_values = etaValues;
stats.errors = errors;
fprintf('test_entropy_prox_eta_zero_limit: %.3e -> %.3e, eta=0 bitwise %d\n', ...
    errors(1), errors(end), stats.eta_zero_bitwise_equal);
end
