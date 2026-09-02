function diagnostic = InteriorKKTResidual(rho, problem, options, gradient)
%INTERIORKKTRESIDUAL Full-target equality KKT and active-set diagnostics.

rho = rho(:);
if nargin < 4 || isempty(gradient)
    gradient = src.discretization.ps.Gradient(rho, problem);
else
    gradient = gradient(:);
end
if nargin < 3 || isempty(options)
    options = struct();
end
if ~isfield(options, 'active_step') || isempty(options.active_step)
    options.active_step = 1;
end
if ~isfield(options, 'active_tol') || isempty(options.active_tol)
    options.active_tol = 1e-12;
end
if ~isfield(options, 'projection_tol') || isempty(options.projection_tol)
    options.projection_tol = 1e-14;
end
if ~isfield(options, 'residual_step') || isempty(options.residual_step)
    options.residual_step = 1;
end

h = problem.grid.h;
onesVector = ones(size(rho));
lambda = -(h * sum(gradient)) / (h * sum(onesVector));
stationarity = gradient + lambda * onesVector;
interiorResidual = sqrt(h * sum(stationarity .^ 2));
[fullPgResidual, fullPgMapping, projectedGradientState] = ...
    src.solvers.FullGradientMapping(rho, gradient, problem, options);
projectedActiveState = src.constraints.ProjectPositiveConservative( ...
    rho - options.active_step * gradient, problem.mass, h, ...
    options.projection_tol);
free = projectedActiveState > options.active_tol;
active = ~free;

diagnostic.gradient = gradient;
diagnostic.lambda = lambda;
diagnostic.stationarity = stationarity;
diagnostic.interior_residual = interiorResidual;
diagnostic.full_pg_residual = fullPgResidual;
diagnostic.full_pg_mapping = fullPgMapping;
diagnostic.projected_gradient_state = projectedGradientState;
diagnostic.projected_active_state = projectedActiveState;
diagnostic.free = free;
diagnostic.active = active;
diagnostic.free_count = nnz(free);
diagnostic.active_count = nnz(active);
diagnostic.projected_free_count = nnz(free);
diagnostic.projected_active_count = nnz(active);
diagnostic.exact_zero_count = nnz(rho == 0);
diagnostic.strict_positive_count = nnz(rho > 0);
diagnostic.mass_residual = src.constraints.Mass(rho, h) - problem.mass;
diagnostic.mass_error = abs(diagnostic.mass_residual);
diagnostic.min_density = min(rho);
end
