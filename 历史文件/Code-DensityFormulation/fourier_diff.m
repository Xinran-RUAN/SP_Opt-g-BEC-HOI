function [dX] = fourier_diff(X, L)
    N = length(X);
    Lambda = 2 * pi * 1i / L * [0:(N/2-1), (-N/2:-1)];
    Lambda = reshape(Lambda, size(X)); % Lambda与X形状保持一致
%%
    FX = fft(X);
    FdX = Lambda .* FX;
    dX = real(ifft(FdX));
end