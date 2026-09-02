function stats = test_smooth_step_cinf()
%TEST_SMOOTH_STEP_CINF Flat endpoints, symmetry, and monotonicity.

t = [-2; 0; linspace(1e-6, 1 - 1e-6, 1001).'; 1; 3];
chi = model.SmoothStepCInf(t);
interior = chi(3:end-2);
stats.endpoint_error = max(abs([chi(1); chi(2); ...
    chi(end-1) - 1; chi(end) - 1]));
stats.symmetry_error = max(abs(interior + flipud(interior) - 1));
stats.minimum_increment = min(diff(chi));
assert(stats.endpoint_error == 0, ...
    'C-infinity step did not preserve exact endpoint plateaus.');
assert(stats.symmetry_error <= 10 * eps, ...
    'C-infinity step symmetry was inconsistent.');
assert(stats.minimum_increment >= -10 * eps, ...
    'C-infinity step must be nondecreasing.');
assert(all(isfinite(chi)) && all(chi >= 0) && all(chi <= 1));
fprintf(['test_smooth_step_cinf: endpoint %.3e, symmetry %.3e, ' ...
    'min increment %.3e\n'], stats.endpoint_error, ...
    stats.symmetry_error, stats.minimum_increment);
end
