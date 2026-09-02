function grid = SetupGrid2D(parameters)
%SETUPGRID2D Periodic tensor grid on [-L,L) x [-L,L).

required = {'L', 'Nx', 'Ny'};
if ~isstruct(parameters) || ~all(isfield(parameters, required))
    error('model:SetupGrid2D:InvalidInput', ...
        'parameters must contain L, Nx, and Ny.');
end
L = parameters.L;
Nx = parameters.Nx;
Ny = parameters.Ny;
if ~isscalar(L) || ~isfinite(L) || L <= 0
    error('model:SetupGrid2D:InvalidL', 'L must be positive and finite.');
end
if ~validCount(Nx) || ~validCount(Ny)
    error('model:SetupGrid2D:InvalidN', ...
        'Nx and Ny must be even integers not smaller than four.');
end

hx = 2 * L / Nx;
hy = 2 * L / Ny;
x = -L + (0:Nx-1)' * hx;
y = -L + (0:Ny-1)' * hy;
[X, Y] = meshgrid(x, y);

grid.dimension = 2;
grid.L = L;
grid.Lx = L;
grid.Ly = L;
grid.Nx = Nx;
grid.Ny = Ny;
grid.N = Nx * Ny;
grid.hx = hx;
grid.hy = hy;
grid.h = hx * hy;
grid.cell_measure = grid.h;
grid.domain_measure = (2 * L) ^ 2;
% Compatibility name used by the dimension-independent solver validation.
grid.domain_length = grid.domain_measure;
grid.shape = [Ny, Nx];
grid.x = x;
grid.y = y;
grid.X = X;
grid.Y = Y;
end

function tf = validCount(N)
tf = isscalar(N) && isfinite(N) && N == round(N) ...
    && N >= 4 && mod(N, 2) == 0;
end
