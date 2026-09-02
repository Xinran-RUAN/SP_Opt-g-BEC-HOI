function value = Value(rho, grid)
%VALUE Discrete entropy H_h(rho) with the convention 0*log(0)=0.

rho = rho(:);
if any(~isfinite(rho)) || ~isreal(rho) || any(rho < 0)
    error('src:entropy:Value:InvalidDensity', ...
        'rho must be a finite, real, nonnegative vector.');
end
if isstruct(grid)
    h = grid.h;
else
    h = grid;
end
positive = rho > 0;
summand = -rho;
summand(positive) = rho(positive) .* (log(rho(positive)) - 1);
value = h * sum(summand);
end
