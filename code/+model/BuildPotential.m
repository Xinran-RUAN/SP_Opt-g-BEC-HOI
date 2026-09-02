function V = BuildPotential(grid, potentialLabel)
%BUILDPOTENTIAL Evaluate the requested potential at collocation points.

if nargin < 2 || isempty(potentialLabel)
    potentialLabel = 'harmonic';
end
switch lower(char(potentialLabel))
    case 'harmonic'
        V = 0.5 * grid.x .^ 2;
    otherwise
        error('model:BuildPotential:UnsupportedPotential', ...
            'Unsupported potential "%s".', potentialLabel);
end
end
