function values = EvaluateFourierSeries1D(coefficients, sourceL, xQuery)
%EVALUATEFOURIERSERIES1D Evaluate the project Fourier series at xQuery.
%
% coefficients use FourierCoefficients' shifted convention
%   c_k = fft(u)_k/N,  k=-N/2,...,N/2-1,
% on the source grid x_j=-L+2Lj/N. Consequently the basis phase is
% exp(i*pi*k*(x+L)/L), including the source-grid origin shift.

coefficients = coefficients(:);
N = numel(coefficients);
if isempty(coefficients) || mod(N, 2) ~= 0 ...
        || any(~isfinite(coefficients))
    error('src:discretization:ps:EvaluateFourierSeries1D:InvalidCoefficients', ...
        'coefficients must be a finite, nonempty, even-length vector.');
end
if ~isscalar(sourceL) || ~isfinite(sourceL) || sourceL <= 0
    error('src:discretization:ps:EvaluateFourierSeries1D:InvalidL', ...
        'sourceL must be a positive finite scalar.');
end
if ~isreal(xQuery) || any(~isfinite(xQuery(:)))
    error('src:discretization:ps:EvaluateFourierSeries1D:InvalidQuery', ...
        'xQuery must contain finite real coordinates.');
end

queryShape = size(xQuery);
x = xQuery(:);
modes = (-N/2:N/2-1).';
values = complex(zeros(size(x)));
% Bound the temporary phase matrix for larger arbitrary query sets.
blockSize = max(1, floor(2e6 / N));
for first = 1:blockSize:numel(x)
    last = min(numel(x), first + blockSize - 1);
    phase = (pi / sourceL) * (x(first:last) + sourceL) * modes.';
    values(first:last) = exp(1i * phase) * coefficients;
end
imaginaryScale = max(1, max(abs(values)));
if max(abs(imag(values))) <= 100 * eps(imaginaryScale)
    values = real(values);
end
values = reshape(values, queryShape);
end
