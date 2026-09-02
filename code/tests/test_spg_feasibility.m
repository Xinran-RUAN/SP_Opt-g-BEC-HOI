function stats = test_spg_feasibility()
%TEST_SPG_FEASIBILITY Check every accepted line-segment iterate.

config = smallSPGConfig();
rng(41);
grid = model.SetupGrid1D(config.parameters);
rho0 = 0.2 + rand(grid.N, 1);
rho0 = rho0 / src.constraints.Mass(rho0, grid.h);
result = experiments.SolveGroundState(config, rho0);
history = result.history.main;

maxMassError = max(history.mass_error);
minDensity = min(history.min_density);
assert(~result.diagnostics.failed, 'SPG reported a line-search failure.');
assert(minDensity >= -1e-13, ...
    'SPG minimum accepted density is %.3e.', minDensity);
assert(maxMassError <= 1e-12, ...
    'SPG maximum accepted mass error is %.3e.', maxMassError);

stats.max_mass_error = maxMassError;
stats.min_density = minDensity;
stats.iterations = result.diagnostics.main_iterations;
fprintf('test_spg_feasibility: min %.3e, mass %.3e\n', ...
    minDensity, maxMassError);
end

function config = smallSPGConfig()
config = experiments.DefaultConfig();
config.parameters.L = 8;
config.parameters.N = 64;
config.solver.name = 'spg';
config.solver.polish_mode = 'none';
config.solver.pg_tol = 1e-30;
config.solver.max_iter = 40;
config.solver.display = false;
end
