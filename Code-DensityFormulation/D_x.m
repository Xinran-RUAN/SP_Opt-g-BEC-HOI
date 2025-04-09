% Given vector x, compute D*x (numerical approx of grad(x)) 
%   where D = 1 / N * Fc * diag(mu) * conj(Fc')
function[Dx] = D_x(x, mu)
Fx = myfft(x);
DFX = 1i * mu .* Fx;
Dx = myifft(DFX);
