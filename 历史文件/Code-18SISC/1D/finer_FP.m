% interpolation via FP
function[Rho_f] = finer_FP(Rho,NN)
for jj = 1 : NN 
    N = length(Rho);
    my_xi = exp(1i*pi/N);
    my_lmd_shift = [ my_xi.^(0:N/2-1), my_xi.^(-N/2:-1)].';
    Rho_shift = real(ifft(my_lmd_shift.*fft(Rho)));
    Rho_mat = [Rho' ; Rho_shift'];
    Rho = Rho_mat(:);
end
Rho_f=Rho;

% Be careful when doing the transpose:
% Two types: ' and .'(will not do the conjugate)