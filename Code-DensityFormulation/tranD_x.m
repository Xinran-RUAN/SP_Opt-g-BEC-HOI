% Given column vector x, compute tran(D)*x 
%   where D = i / N * Fc * diag(mu) * conj(tran(Fc))
%   tran(D) = i / N * conj(Fc) * diag(mu) * tran(Fc)
function[tranDx] = tranD_x(x, mu)
n = length(x);
% 1 / N * tran(Fc)
f1 = ifft(x);
tFx_n = [f1(n/2+1:n); f1(1:n/2)];
% 1i * diag  
DFX = 1i * mu .* tFx_n;
% conj(Fc)
DFX_p = [DFX(n/2+1:n); DFX(1:n/2)];
tranDx = fft(DFX_p);
