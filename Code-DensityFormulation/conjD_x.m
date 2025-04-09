% Given vector x, compute conj(D)*x, 
%   where D = 1 / N * Fc * diag(mu) * conj(Fc')
function[conjDx] = conjD_x(x, mu)
Fx = myfft(x);
DFX = 1i * mu .* Fx;
conjDx = conj(myifft(DFX));