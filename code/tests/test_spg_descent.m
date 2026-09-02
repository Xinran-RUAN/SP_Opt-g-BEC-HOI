function stats = test_spg_descent()
%TEST_SPG_DESCENT Verify h-inner-product descent of projected directions.

config = experiments.DefaultConfig();
config.parameters.L = 8;
config.parameters.N = 64;
config.solver.name = 'spg';
config.solver.polish_mode = 'none';
config.solver.pg_tol = 1e-30;
config.solver.max_iter = 40;
config.solver.display = false;
result = experiments.SolveGroundState(config);
descent = result.history.main.descent;
maximumDescent = max(descent);
assert(maximumDescent <= 1e-12, ...
    'SPG projected direction has positive descent diagnostic %.3e.', ...
    maximumDescent);

stats.maximum_descent = maximumDescent;
stats.minimum_descent = min(descent);
fprintf('test_spg_descent: maximum <G,d>_h %.3e\n', maximumDescent);
end
