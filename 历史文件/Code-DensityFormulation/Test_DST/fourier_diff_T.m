function [drho_T] = fourier_diff_T(rho, Lambda)
    iFrho = ifft(rho);
    iFdrho = Lambda .* iFrho;
    drho_T = real(fft(iFdrho));
end