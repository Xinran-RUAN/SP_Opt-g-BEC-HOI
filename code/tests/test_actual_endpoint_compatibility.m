function stats = test_actual_endpoint_compatibility()
%TEST_ACTUAL_ENDPOINT_COMPATIBILITY Analytic harmonic jump on FFT state.

epsilon = 1e-3;
sigma = 1e-6;
L = 8;
parameters.L = L;
parameters.N = 128;
grid = model.SetupGrid1D(parameters);
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = 0.5 * grid.x .^ 2;
problem.beta = 10;
problem.delta = 10;
problem.mass = 1;
problem.fisher_regularization.epsilon = epsilon;
problem.fisher_regularization.s_epsilon = @(rho) rho + epsilon;
problem.fisher_regularization.ds_epsilon = @(rho) ones(size(rho));
problem.fisher_regularization.d2s_epsilon = @(rho) zeros(size(rho));
problem.potential_regularization.sigma = sigma;
problem.potential_regularization.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
problem.potential_regularization.dp_sigma = @(rho) rho ./ hypot(rho, sigma);
problem.potential_regularization.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
rho = 0.1 + 0.01 * cos(pi * grid.x / L);
endpoint.V = [0.5 * L ^ 2, 0.5 * L ^ 2];
endpoint.dV = [-L, L];
diagnostic = src.diagnostics.ActualEndpointCompatibility( ...
    rho, problem, endpoint);
expected = 2 * L * problem.potential_regularization.dp_sigma(rho(1));
stats.rho_value_jump = diagnostic.rho_value_jump;
stats.rho_derivative_jump = diagnostic.rho_derivative_jump;
stats.GV_derivative_error = abs(diagnostic.GV_derivative_jump - expected);
assert(stats.rho_value_jump == 0);
assert(stats.rho_derivative_jump == 0);
assert(stats.GV_derivative_error <= 1e-13);
fprintf('test_actual_endpoint_compatibility: GV jump error %.3e\n', ...
    stats.GV_derivative_error);
end
