function stats = test_fourier_derivative()
%TEST_FOURIER_DERIVATIVE Spectral accuracy and skew-adjointness of D.

parameters = model.DefaultParameters1D();
parameters.L = 7;
parameters.N = 64;
grid = model.SetupGrid1D(parameters);
plan = src.discretization.ps.Plan1D(grid);
x = grid.x;

u = sin(3 * pi * x / grid.L) + 0.4 * cos(5 * pi * x / grid.L);
duExact = (3 * pi / grid.L) * cos(3 * pi * x / grid.L) ...
    - 0.4 * (5 * pi / grid.L) * sin(5 * pi * x / grid.L);
du = src.discretization.ps.FirstDerivative(u, plan);

rng(11);
v = randn(grid.N, 1);
dv = src.discretization.ps.FirstDerivative(v, plan);
adjointDefect = abs(grid.h * sum(du .* v + u .* dv));
derivativeError = max(abs(du - duExact));

assert(derivativeError <= 2e-12, ...
    'Fourier derivative error %.3e is too large.', derivativeError);
assert(adjointDefect <= 2e-12, ...
    'Skew-adjoint defect %.3e is too large.', adjointDefect);

stats.derivative_error = derivativeError;
stats.adjoint_defect = adjointDefect;
fprintf('test_fourier_derivative: derivative %.3e, adjoint %.3e\n', ...
    derivativeError, adjointDefect);
end
