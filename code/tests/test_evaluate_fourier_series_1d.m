function stats = test_evaluate_fourier_series_1d()
%TEST_EVALUATE_FOURIER_SERIES_1D Nodal reconstruction and arbitrary points.

L = 7;
N = 128;
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
values = 1 + 0.2 * cos(3 * pi * grid.x / L) ...
    - 0.07 * sin(11 * pi * grid.x / L);
[coefficients, ~] = src.discretization.ps.FourierCoefficients(values);
nodal = src.discretization.ps.EvaluateFourierSeries1D( ...
    coefficients, L, grid.x);
xQuery = linspace(-5, 5, 1001).';
exact = 1 + 0.2 * cos(3 * pi * xQuery / L) ...
    - 0.07 * sin(11 * pi * xQuery / L);
arbitrary = src.discretization.ps.EvaluateFourierSeries1D( ...
    coefficients, L, xQuery);
stats.nodal_error = max(abs(nodal - values));
stats.arbitrary_error = max(abs(arbitrary - exact));
assert(stats.nodal_error <= 1e-12, ...
    'Fourier arbitrary-point evaluator failed nodal reconstruction.');
assert(stats.arbitrary_error <= 1e-12, ...
    'Fourier arbitrary-point evaluator used inconsistent phases.');
fprintf('test_evaluate_fourier_series_1d: nodal %.3e, arbitrary %.3e\n', ...
    stats.nodal_error, stats.arbitrary_error);
end
