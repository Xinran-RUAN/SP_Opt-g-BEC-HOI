function stats = test_transition_layer_diagnostics()
%TEST_TRANSITION_LAYER_DIAGNOSTICS Crossing interpolation and symmetry.

L = 8;
N = 1024;
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
plan = src.discretization.ps.Plan1D(grid);
q = 0.05 + 0.45 * (1 + cos(pi * grid.x / L));
sigma = 1e-3;
rho = sigma * q ./ sqrt(1 - q .^ 2);
potential.sigma = sigma;
potential.dp_sigma = @(r) r ./ hypot(r, sigma);
potential.d2p_sigma = @(r) ...
    (sigma ./ hypot(r, sigma)) .^ 2 ./ hypot(r, sigma);
diagnostic = src.diagnostics.TransitionLayerDiagnostics( ...
    rho, grid, plan, potential);

stats.valid = diagnostic.primary.valid;
stats.asymmetry = diagnostic.primary.asymmetry;
stats.width = diagnostic.primary.width_mean;
stats.points = diagnostic.primary.points_per_layer;
assert(stats.valid, 'Manufactured q layer crossings were not found.');
assert(stats.asymmetry <= 1e-12);
assert(stats.width > 0 && stats.points > 0);
fprintf(['test_transition_layer_diagnostics: width %.3e, points %.2f, ' ...
    'asymmetry %.3e\n'], stats.width, stats.points, stats.asymmetry);
end
