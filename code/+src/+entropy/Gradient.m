function gradient = Gradient(rho)
%GRADIENT Entropy gradient under <u,v>_h = h*sum(u.*v).

rho = rho(:);
if any(~isfinite(rho)) || ~isreal(rho) || any(rho < 0)
    error('src:entropy:Gradient:InvalidDensity', ...
        'rho must be a finite, real, nonnegative vector.');
end
gradient = log(rho);
end
