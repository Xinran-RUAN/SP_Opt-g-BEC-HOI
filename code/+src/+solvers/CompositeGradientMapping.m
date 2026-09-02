function [residual, mapping, proximalState, proxInfo] = ...
    CompositeGradientMapping(rho, problem, solver)
%COMPOSITEGRADIENTMAPPING Fixed-step potential-prox stationarity mapping.

if entropyActive(problem)
    error('src:solvers:CompositeGradientMapping:EntropyConflict', ...
        'Potential composite mapping is separate from the entropy branch.');
end
tau = solver.residual_step;
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('src:solvers:CompositeGradientMapping:InvalidStep', ...
        'solver.residual_step must be positive and finite.');
end
smoothGradient = src.discretization.ps.SmoothGradient(rho, problem);
proxSolver = solver;
proxSolver.splitting = 'potential_prox';
[proximalState, proxInfo] = src.solvers.ApplyProx( ...
    rho(:) - tau * smoothGradient(:), tau, problem, proxSolver);
mapping = (rho(:) - proximalState(:)) / tau;
residual = sqrt(problem.grid.h * sum(mapping .^ 2));
end

function active = entropyActive(problem)
active = isfield(problem, 'entropy') ...
    && isfield(problem.entropy, 'enabled') && problem.entropy.enabled ...
    && isfield(problem.entropy, 'eta') && problem.entropy.eta > 0;
end
