function mass = Mass(rho, h)
%MASS Discrete collocation mass h*sum(rho).

mass = h * sum(rho(:));
end
