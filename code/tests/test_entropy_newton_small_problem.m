function stats = test_entropy_newton_small_problem()
%TEST_ENTROPY_NEWTON_SMALL_PROBLEM Verify interior KKT improvement.

config = experiments.DefaultConfig();
config.parameters.L = 2;
config.parameters.N = 32;
config.entropy.enabled = true;
config.entropy.eta = 0.1;
config.solver.name = 'fista_cd';
config.solver.pg_tol = 1e-6;
config.solver.max_iter = 10000;
config.solver.polish_mode = 'always';
config.solver.polish.entry_pg_tol = 1e-3;
config.solver.polish.pg_tol = 1e-8;
config.solver.display = false;
result = experiments.SolveGroundState(config);
before = result.diagnostics.kkt_before_polish;
after = result.diagnostics.kkt_after_polish;
assert(result.diagnostics.polish_attempted, ...
    'Entropy Newton was not attempted.');
assert(~result.diagnostics.polish_failed, ...
    'Entropy Newton failed: %s', result.diagnostics.polish_failure_message);
assert(after <= min(1e-8, 0.01 * before), ...
    'Entropy Newton did not reduce KKT residual: %.3e -> %.3e.', ...
    before, after);
assert(result.diagnostics.mass_error <= 1e-12 && min(result.rho) > 0, ...
    'Entropy Newton lost feasibility/interiority.');
stats.kkt_before = before;
stats.kkt_after = after;
stats.pg_after = result.diagnostics.pg_after_polish;
stats.iterations = result.diagnostics.polish_iterations;
fprintf('test_entropy_newton_small_problem: KKT %.3e -> %.3e\n', ...
    before, after);
end
