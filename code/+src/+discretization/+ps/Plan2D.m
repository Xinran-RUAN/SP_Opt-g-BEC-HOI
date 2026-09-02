function plan = Plan2D(grid)
%PLAN2D FFT differentiation plan for a periodic even tensor grid.
%
% Arrays use MATLAB's Ny-by-Nx orientation: rows correspond to y and
% columns to x.  Nyquist derivative multipliers are set to zero in both
% directions, matching the production 1D skew-adjoint convention.

required = {'Nx', 'Ny', 'Lx', 'Ly', 'hx', 'hy', 'N'};
if ~isstruct(grid) || ~all(isfield(grid, required))
    error('src:discretization:ps:Plan2D:InvalidGrid', ...
        'grid must contain Nx, Ny, Lx, Ly, hx, hy, and N.');
end
if grid.N ~= grid.Nx * grid.Ny ...
        || mod(grid.Nx, 2) ~= 0 || mod(grid.Ny, 2) ~= 0
    error('src:discretization:ps:Plan2D:InvalidSize', ...
        'Require grid.N=Nx*Ny with even Nx and Ny.');
end

integerModesX = [0:(grid.Nx/2-1), 0, (-grid.Nx/2+1):-1];
integerModesY = [0:(grid.Ny/2-1), 0, (-grid.Ny/2+1):-1]';
kx = (pi / grid.Lx) * integerModesX;
ky = (pi / grid.Ly) * integerModesY;

plan.dimension = 2;
plan.N = grid.N;
plan.Nx = grid.Nx;
plan.Ny = grid.Ny;
plan.Lx = grid.Lx;
plan.Ly = grid.Ly;
plan.hx = grid.hx;
plan.hy = grid.hy;
plan.h = grid.h;
plan.shape = [grid.Ny, grid.Nx];
plan.integer_modes_x = integerModesX;
plan.integer_modes_y = integerModesY;
plan.wavenumbers_x = kx;
plan.wavenumbers_y = ky;
plan.first_derivative_multiplier_x = 1i * kx;
plan.first_derivative_multiplier_y = 1i * ky;
plan.coefficient_convention = 'c_k = fft2(u)_k / (Nx*Ny)';
plan.nyquist_derivative_is_zero = true;
end
