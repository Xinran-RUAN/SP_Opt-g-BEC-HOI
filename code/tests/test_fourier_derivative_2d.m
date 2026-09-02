function stats = test_fourier_derivative_2d()
%TEST_FOURIER_DERIVATIVE_2D Tensor FFT derivatives and adjoint identity.

parameters.L = 3;
parameters.Nx = 32;
parameters.Ny = 40;
grid = model.SetupGrid2D(parameters);
plan = src.discretization.ps.Plan2D(grid);
u = sin(pi * grid.X / grid.L) ...
    + 0.3 * cos(3 * pi * grid.Y / grid.L) ...
    + 0.1 * sin(2 * pi * grid.X / grid.L) ...
    .* cos(4 * pi * grid.Y / grid.L);
ux = (pi / grid.L) * cos(pi * grid.X / grid.L) ...
    + 0.2 * (pi / grid.L) * cos(2 * pi * grid.X / grid.L) ...
    .* cos(4 * pi * grid.Y / grid.L);
uy = -0.9 * (pi / grid.L) * sin(3 * pi * grid.Y / grid.L) ...
    - 0.4 * (pi / grid.L) * sin(2 * pi * grid.X / grid.L) ...
    .* sin(4 * pi * grid.Y / grid.L);
computed = src.discretization.ps.SpatialGradient(u(:), plan);
stats.x_error = max(abs(computed(:, 1) - ux(:)));
stats.y_error = max(abs(computed(:, 2) - uy(:)));

rng(7);
a = randn(grid.N, 1);
b = randn(grid.N, 2);
left = grid.h * sum(src.discretization.ps.SpatialGradient(a, plan) .* b, 'all');
right = grid.h * sum(a .* src.discretization.ps.AdjointSpatialGradient(b, plan));
stats.adjoint_error = abs(left - right) / max([1, abs(left), abs(right)]);
assert(max(stats.x_error, stats.y_error) <= 1e-12, ...
    '2D derivative error is %.3e.', max(stats.x_error, stats.y_error));
assert(stats.adjoint_error <= 1e-12, ...
    '2D derivative adjoint error is %.3e.', stats.adjoint_error);
fprintf('test_fourier_derivative_2d: Dx %.3e, Dy %.3e, adjoint %.3e\n', ...
    stats.x_error, stats.y_error, stats.adjoint_error);
end
