function state = EvaluateState( ...
    rho, problem, solver, energy, gradient, computeProjectedResidual)
%EVALUATESTATE Objective, unified stationarity, and feasibility diagnostics.

rho = rho(:);
if nargin < 4 || isempty(energy)
    [energy, state.energy_parts] = ...
        src.discretization.ps.Energy(rho, problem);
else
    state.energy_parts = struct();
end
if nargin < 5 || isempty(gradient)
    gradient = src.discretization.ps.Gradient(rho, problem);
end
if nargin < 6 || isempty(computeProjectedResidual)
    computeProjectedResidual = true;
end
state.energy = energy;
state.physical_energy = energy;
if entropyActive(problem)
    state.entropy_value = src.entropy.Value(rho, problem.grid);
    state.augmented_energy = energy ...
        + problem.entropy.eta * state.entropy_value;
else
    state.entropy_value = src.entropy.Value(rho, problem.grid);
    state.augmented_energy = energy;
end
state.gradient = gradient;
if computeProjectedResidual
    [state.pg_residual, state.pg_mapping, state.projected_state] = ...
        src.solvers.ProjectedGradientResidual(rho, gradient, problem, solver);
else
    state.pg_residual = NaN;
    state.pg_mapping = [];
    state.projected_state = [];
end
if entropyActive(problem)
    state.kkt = src.entropy.KKTResidual(rho, gradient, problem);
else
    state.kkt = src.solvers.KKTResidual( ...
        rho, gradient, problem, solver.active_tol);
end
state.kkt_residual = state.kkt.kkt_residual;
state.mass_error = abs(src.constraints.Mass(rho, problem.grid.h) - problem.mass);
state.min_density = min(rho);
end

function active = entropyActive(problem)
active = isfield(problem, 'entropy') ...
    && isfield(problem.entropy, 'enabled') && problem.entropy.enabled ...
    && isfield(problem.entropy, 'eta') && problem.entropy.eta > 0;
end
