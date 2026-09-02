function [dX_T] = fourier_diff_T(X, L)
%%
    N = length(X);
    Lambda = 2 * pi * 1i / L * [0:(N/2-1), (-N/2:-1)];
    Lambda = reshape(Lambda, size(X)); % Lambda与X形状保持一致
%%
    iFX = ifft(X);
    iFdX = Lambda .* iFX;
    dX_T = real(fft(iFdX));
end