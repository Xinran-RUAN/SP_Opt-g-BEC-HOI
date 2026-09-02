function stats = test_pdas_small_problem()
%TEST_PDAS_SMALL_PROBLEM Verify local stationarity improvement and feasibility.

config = experiments.DefaultConfig();
config.parameters.N = 32;
config.solver.name = 'spg';
config.solver.pg_tol = 1e-8;
config.solver.polish_mode = 'always';
config.solver.polish.pg_tol = 1e-12;
config.solver.display = false;
result = experiments.SolveGroundState(config);
before = result.diagnostics.pg_before_polish;
after = result.diagnostics.pg_after_polish;
energyChange = result.diagnostics.polish_energy_change;

assert(result.diagnostics.polish_attempted, 'PDAS polish was not attempted.');
assert(result.diagnostics.polish_converged, 'PDAS polish did not converge.');
assert(after <= min(1e-11, 0.01 * before), ...
    'PDAS did not substantially reduce PG residual: %.3e -> %.3e.', ...
    before, after);
assert(energyChange <= 1e-10, ...
    'PDAS materially increased energy by %.3e.', energyChange);
assert(result.diagnostics.mass_error <= 1e-12, ...
    'PDAS mass error is %.3e.', result.diagnostics.mass_error);
assert(result.diagnostics.min_density >= -1e-13, ...
    'PDAS minimum density is %.3e.', result.diagnostics.min_density);

stats.pg_before = before;
stats.pg_after = after;
stats.kkt_after = result.diagnostics.kkt_after_polish;
stats.energy_change = energyChange;
stats.iterations = result.diagnostics.polish_iterations;
fprintf('test_pdas_small_problem: PG %.3e -> %.3e, KKT %.3e\n', ...
    before, after, stats.kkt_after);
end
