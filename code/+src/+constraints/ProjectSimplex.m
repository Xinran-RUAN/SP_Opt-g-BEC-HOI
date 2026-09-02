function rho = ProjectSimplex(z, mass, h)
%PROJECTSIMPLEX Euclidean projection onto rho>=0, h*sum(rho)=mass.

if ~isscalar(mass) || ~isfinite(mass) || mass <= 0 || ...
        ~isscalar(h) || ~isfinite(h) || h <= 0
    error('src:constraints:ProjectSimplex:InvalidParameters', ...
        'mass and h must be positive finite scalars.');
end
wasRow = isrow(z);
z = z(:);
if isempty(z) || any(~isfinite(z)) || ~isreal(z)
    error('src:constraints:ProjectSimplex:InvalidInput', ...
        'z must be a nonempty finite real vector.');
end

targetSum = mass / h;
u = sort(z, 'descend');
cssv = cumsum(u) - targetSum;
index = (1:numel(z))';
active = find(u - cssv ./ index > 0, 1, 'last');
if isempty(active)
    error('src:constraints:ProjectSimplex:ProjectionFailure', ...
        'Could not determine the active simplex set.');
end
threshold = cssv(active) / active;
rho = max(z - threshold, 0);
if wasRow
    rho = rho.';
end
end
