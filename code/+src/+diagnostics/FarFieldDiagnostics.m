function diagnostic = FarFieldDiagnostics(rho, grid, plan)
%FARFIELDDIAGNOSTICS Density and first-derivative behavior near box edges.

rho = rho(:);
if nargin < 3 || isempty(plan)
    plan = src.discretization.ps.Plan1D(grid);
end
if numel(rho) ~= grid.N
    error('src:diagnostics:FarFieldDiagnostics:SizeMismatch', ...
        'rho must contain grid.N entries.');
end
tailWidth = min(2, 0.1 * grid.L);
tailMask = abs(grid.x) >= grid.L - tailWidth;
boundaryStrip = abs(grid.x) >= 0.8 * grid.L;
derivative = src.discretization.ps.FirstDerivative(rho, plan);

diagnostic.tail_width = tailWidth;
diagnostic.tail_mass = grid.h * sum(rho(tailMask));
diagnostic.tail_max = max(rho(tailMask));
diagnostic.boundary_strip_mass = grid.h * sum(rho(boundaryStrip));
diagnostic.boundary_strip_max = max(rho(boundaryStrip));
diagnostic.min_density = min(rho);
diagnostic.edge_density_left = rho(1);
diagnostic.edge_density_right = rho(end);
diagnostic.edge_density = max(abs([rho(1), rho(end)]));
diagnostic.edge_abs_drho_left = abs(derivative(1));
diagnostic.edge_abs_drho_right = abs(derivative(end));
diagnostic.edge_abs_drho = max(abs([derivative(1), derivative(end)]));
diagnostic.boundary_strip_max_abs_drho = ...
    max(abs(derivative(boundaryStrip)));
end
