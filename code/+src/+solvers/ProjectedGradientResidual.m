function [residual, mapping, projected] = ProjectedGradientResidual( ...
    rho, gradient, problem, solver)
%PROJECTEDGRADIENTRESIDUAL Unified fixed-step stationarity residual.
%
% Every first-order solver and the PDAS polish use the same residual step,
% tau_r = solver.residual_step, so their reported residuals are comparable.

tau = solver.residual_step;
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('src:solvers:ProjectedGradientResidual:InvalidStep', ...
        'solver.residual_step must be positive and finite.');
end
if ~entropyActive(problem) && isfield(solver, 'splitting') ...
        && strcmpi(solver.splitting, 'potential_prox')
    [residual, mapping, projected] = ...
        src.solvers.CompositeGradientMapping(rho, problem, solver);
    return;
end
projected = src.solvers.ApplyProx( ...
    rho(:) - tau * gradient(:), tau, problem, solver);
mapping = (rho(:) - projected(:)) / tau;
residual = sqrt(problem.grid.h * sum(mapping .^ 2));
end

function active = entropyActive(problem)
active = isfield(problem, 'entropy') ...
    && isfield(problem.entropy, 'enabled') && problem.entropy.enabled ...
    && isfield(problem.entropy, 'eta') && problem.entropy.eta > 0;
end
