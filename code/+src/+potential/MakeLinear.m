function potential = MakeLinear()
%MAKELINEAR Paper baseline p_0(rho)=rho in the generic handle interface.

potential.sigma = 0;
potential.p_sigma = @(rho) rho;
potential.dp_sigma = @(rho) ones(size(rho));
potential.d2p_sigma = @(rho) zeros(size(rho));
potential.label = 'p_0(rho) = rho';
potential.prox_type = 'linear';
potential.name = 'linear';
potential.power = 0;
end
