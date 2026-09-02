function [coefficients, modes] = FourierCoefficients(u)
%FOURIERCOEFFICIENTS Return fftshifted c_k=fft(u)_k/N and integer modes.

u = u(:);
N = numel(u);
if mod(N, 2) ~= 0
    error('src:discretization:ps:FourierCoefficients:InvalidN', ...
        'The number of samples must be even.');
end
coefficients = fftshift(fft(u) / N);
modes = (-N/2:N/2-1)';
end
