function [rho, info] = ProxStep(y, gradient, L, problem, solver)
%PROXSTEP Apply the configured composite prox to y-gradient/L.

if ~isscalar(L) || ~isfinite(L) || L <= 0
    error('src:solvers:ProxStep:InvalidL', 'L must be positive and finite.');
end
trial = y(:) - gradient(:) / L;
[rho, info] = src.solvers.ApplyProx(trial, 1 / L, problem, solver);
end
