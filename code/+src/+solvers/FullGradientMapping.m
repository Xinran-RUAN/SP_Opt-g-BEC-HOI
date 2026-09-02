function [residual, mapping, projected] = ...
    FullGradientMapping(rho, fullGradient, problem, solver)
%FULLGRADIENTMAPPING Bo Lin projection of the complete target gradient.

tau = solver.residual_step;
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('src:solvers:FullGradientMapping:InvalidStep', ...
        'solver.residual_step must be positive and finite.');
end
projected = src.constraints.ProjectPositiveConservative( ...
    rho(:) - tau * fullGradient(:), problem.mass, problem.grid.h, ...
    solver.projection_tol);
mapping = (rho(:) - projected(:)) / tau;
residual = sqrt(problem.grid.h * sum(mapping .^ 2));
end
