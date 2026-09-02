function dtu = AdjointFirstDerivative(u, plan)
%ADJOINTFIRSTDERIVATIVE Adjoint of the pseudospectral first derivative.
%
% With the Nyquist derivative set to zero, D^T = -D exactly up to FFT
% roundoff on the periodic uniform grid.

dtu = -src.discretization.ps.FirstDerivative(u, plan);
end
