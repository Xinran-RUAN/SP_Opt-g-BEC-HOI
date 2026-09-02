function du = FirstDerivative(u, plan)
%FIRSTDERIVATIVE Fourier pseudospectral first derivative.

wasRow = isrow(u);
u = u(:);
if numel(u) ~= plan.N
    error('src:discretization:ps:FirstDerivative:SizeMismatch', ...
        'Input length must equal plan.N.');
end
du = real(ifft(plan.first_derivative_multiplier .* fft(u)));
if wasRow
    du = du.';
end
end
