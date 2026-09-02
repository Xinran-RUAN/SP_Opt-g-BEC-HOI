function diagnostic = KKTResidual(rho, gradient, problem, activeTolerance)
%KKTRESIDUAL Direct free/active-set KKT diagnostic at a feasible state.

if nargin < 4 || isempty(activeTolerance)
    activeTolerance = 1e-12;
end
rho = rho(:);
gradient = gradient(:);
free = rho > activeTolerance;
if ~any(free)
    [~, index] = max(rho);
    free(index) = true;
end
active = ~free;
lambda = -mean(gradient(free));
dualStationarity = gradient + lambda;
h = problem.grid.h;

freeResidual = sqrt(h * sum(dualStationarity(free) .^ 2));
dualViolation = min(dualStationarity(active), 0);
dualResidual = sqrt(h * sum(dualViolation .^ 2));
massError = abs(src.constraints.Mass(rho, h) - problem.mass);
negativityVector = min(rho, 0);
negativity = sqrt(h * sum(negativityVector .^ 2));

diagnostic.kkt_residual = sqrt(freeResidual ^ 2 + dualResidual ^ 2 ...
    + massError ^ 2 + negativity ^ 2);
diagnostic.free_residual = freeResidual;
diagnostic.dual_residual = dualResidual;
diagnostic.mass_error = massError;
diagnostic.negativity = negativity;
diagnostic.lambda = lambda;
diagnostic.free = free;
diagnostic.active = active;
diagnostic.free_count = nnz(free);
diagnostic.active_count = nnz(active);
end
