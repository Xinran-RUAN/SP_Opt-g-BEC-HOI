function diagnostic = KKTResidual(rho, gradientPhysical, problem, lambda)
%KKTRESIDUAL Equality-constrained interior KKT diagnostic for eta>0.

rho = rho(:);
gradientPhysical = gradientPhysical(:);
eta = problem.entropy.eta;
massError = abs(src.constraints.Mass(rho, problem.grid.h) - problem.mass);
diagnostic.mass_error = massError;
diagnostic.min_rho = min(rho);
diagnostic.underflow_count = nnz(rho == 0);

if any(rho <= 0)
    diagnostic.min_log_rho = -Inf;
    diagnostic.lambda = NaN;
    diagnostic.stationarity_residual = Inf;
    diagnostic.kkt_residual = Inf;
    diagnostic.status = 'entropy_numerical_underflow';
    return;
end

logRho = log(rho);
if nargin < 4 || isempty(lambda)
    lambda = -mean(gradientPhysical + eta * logRho);
end
stationarity = gradientPhysical + eta * logRho + lambda;
stationarityResidual = sqrt(problem.grid.h * sum(stationarity .^ 2));
diagnostic.min_log_rho = min(logRho);
diagnostic.lambda = lambda;
diagnostic.stationarity_residual = stationarityResidual;
diagnostic.kkt_residual = max(stationarityResidual, massError);
diagnostic.status = 'ok';
end
