function gradient = SpatialGradient(u, plan)
%SPATIALGRADIENT Matrix of Fourier derivatives, one column per direction.

if isfield(plan, 'dimension')
    dimension = plan.dimension;
else
    dimension = 1;
end
gradient = zeros(plan.N, dimension);
for direction = 1:dimension
    gradient(:, direction) = ...
        src.discretization.ps.DirectionalDerivative(u, plan, direction);
end
end
