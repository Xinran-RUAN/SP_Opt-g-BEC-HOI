function rho = InitialDensity(grid, mass, projectionName)
%INITIALDENSITY Gaussian density projected onto the discrete feasible set.

if nargin < 2 || isempty(mass)
    mass = 1;
end
if nargin < 3 || isempty(projectionName)
    projectionName = 'simplex';
end

rho = exp(-grid.x .^ 2) / sqrt(pi);
switch lower(char(projectionName))
    case 'simplex'
        rho = src.constraints.ProjectSimplex(rho, mass, grid.h);
    case 'semismooth'
        rho = src.constraints.ProjectPositiveConservative( ...
            rho, mass, grid.h, 1e-13);
    otherwise
        error('model:InitialDensity:UnknownProjection', ...
            'Unknown projection backend "%s".', projectionName);
end
end
