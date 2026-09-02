function adjoint = AdjointSpatialGradient(field, plan)
%ADJOINTSPATIALGRADIENT Adjoint of the Fourier spatial-gradient operator.

if isfield(plan, 'dimension')
    dimension = plan.dimension;
else
    dimension = 1;
end
if isvector(field) && dimension == 1
    field = field(:);
end
if ~isequal(size(field), [plan.N, dimension])
    error('src:discretization:ps:AdjointSpatialGradient:SizeMismatch', ...
        'field must have size plan.N-by-plan.dimension.');
end
adjoint = zeros(plan.N, 1);
for direction = 1:dimension
    % Each production Fourier derivative is skew-adjoint because its
    % Nyquist multiplier is zero.
    adjoint = adjoint - src.discretization.ps.DirectionalDerivative( ...
        field(:, direction), plan, direction);
end
end
