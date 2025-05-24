function [drho] = fourier_diff(rho, Lambda)
    Frho = fft(rho);
    Fdrho = Lambda .* Frho;
    drho = real(ifft(Fdrho));
end