function values = ValuesFromCoefficients(coefficients)
%VALUESFROMCOEFFICIENTS Nodal values from fftshifted c_k=fft(u)_k/N.

coefficients = coefficients(:);
N = numel(coefficients);
if isempty(coefficients) || mod(N, 2) ~= 0 ...
        || any(~isfinite(coefficients))
    error('src:discretization:ps:ValuesFromCoefficients:InvalidInput', ...
        'coefficients must be a finite, nonempty, even-length vector.');
end
values = ifft(N * ifftshift(coefficients));
imaginaryScale = max(1, max(abs(values)));
if max(abs(imag(values))) <= 100 * eps(imaginaryScale)
    values = real(values);
end
end
