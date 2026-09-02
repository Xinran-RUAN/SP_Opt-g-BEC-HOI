function stats = test_general_potential_prox()
%TEST_GENERAL_POTENTIAL_PROX Generic p_sigma handles, not a named model.

rng(412);
N = 48;
L = 4;
h = 2 * L / N;
x = -L + (0:N-1)' * h;
z = 0.2 * randn(N, 1);
V = 0.1 + x .^ 2;
tau = 0.07;
mass = 1;
potential.sigma = 0.2;
potential.p_sigma = @(rho) 0.5 * rho .^ 2;
potential.dp_sigma = @(rho) rho;
potential.d2p_sigma = @(rho) ones(size(rho));
potential.label = 'p_sigma(rho) = rho^2/2';
potential.prox_type = 'generic_convex';
options.mass_tol = 1e-13;
options.inner_tol = 1e-13;
[rho, info] = src.potential.PositiveConservativeProx( ...
    z, tau, V, potential, mass, h, options);
stats.mass_error = info.mass_error;
stats.kkt_residual = info.kkt_residual;
stats.min_density = min(rho);
assert(stats.mass_error <= 1e-12 && stats.min_density >= -1e-14 ...
    && stats.kkt_residual <= 2e-12);
fprintf('test_general_potential_prox: mass %.3e, KKT %.3e, min %.3e\n', ...
    stats.mass_error, stats.kkt_residual, stats.min_density);
end
