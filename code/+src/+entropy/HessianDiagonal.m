function diagonal = HessianDiagonal(rho)
%HESSIANDIAGONAL Diagonal 1/rho of the entropy Hessian.

rho = rho(:);
if any(~isfinite(rho)) || ~isreal(rho) || any(rho < 0)
    error('src:entropy:HessianDiagonal:InvalidDensity', ...
        'rho must be a finite, real, nonnegative vector.');
end
diagonal = 1 ./ rho;
end
