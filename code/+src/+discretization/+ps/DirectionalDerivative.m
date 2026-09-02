function derivative = DirectionalDerivative(u, plan, direction)
%DIRECTIONALDERIVATIVE Fourier derivative along one spatial direction.

u = u(:);
if numel(u) ~= plan.N
    error('src:discretization:ps:DirectionalDerivative:SizeMismatch', ...
        'Input length must equal plan.N.');
end
dimension = spatialDimension(plan);
if direction ~= round(direction) || direction < 1 || direction > dimension
    error('src:discretization:ps:DirectionalDerivative:Direction', ...
        'direction must be an integer between one and plan.dimension.');
end
if dimension == 1
    derivative = src.discretization.ps.FirstDerivative(u, plan);
    return;
end

U = reshape(u, plan.Ny, plan.Nx);
Uhat = fft2(U);
if direction == 1
    multiplier = plan.first_derivative_multiplier_x;
else
    multiplier = plan.first_derivative_multiplier_y;
end
derivative = real(ifft2(multiplier .* Uhat));
derivative = derivative(:);
end

function dimension = spatialDimension(plan)
if isfield(plan, 'dimension')
    dimension = plan.dimension;
else
    dimension = 1;
end
end
