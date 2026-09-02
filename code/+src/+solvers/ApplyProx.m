function [rho, info] = ApplyProx(z, tau, problem, solver)
%APPLYPROX Composite prox with a bitwise-compatible eta=0 projection path.

if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('src:solvers:ApplyProx:InvalidStep', ...
        'tau must be a positive finite scalar.');
end
entropyActive = isfield(problem, 'entropy') ...
    && isfield(problem.entropy, 'enabled') && problem.entropy.enabled ...
    && isfield(problem.entropy, 'eta') && problem.entropy.eta > 0;
if entropyActive
    [rho, info] = src.entropy.SimplexProx( ...
        z, tau * problem.entropy.eta, problem.mass, ...
        problem.grid.h, problem.entropy.prox);
    info.backend = 'entropy_simplex';
    return;
end

if isfield(solver, 'splitting') ...
        && strcmpi(solver.splitting, 'potential_prox')
    if ~isfield(problem, 'potential_regularization') ...
            || ~isfield(problem, 'V')
        error('src:solvers:ApplyProx:MissingPotentialProblem', ...
            'potential_prox requires problem.V and potential_regularization.');
    end
    if isfield(solver, 'potential_prox')
        options = solver.potential_prox;
    else
        options = struct();
    end
    options.projection_tol = solver.projection_tol;
    [rho, info] = src.potential.PositiveConservativeProx( ...
        z, tau, problem.V, problem.potential_regularization, ...
        problem.mass, problem.grid.h, options);
    return;
end

switch lower(char(solver.projection_name))
    case 'simplex'
        rho = src.constraints.ProjectSimplex(z, problem.mass, problem.grid.h);
    case 'semismooth'
        rho = src.constraints.ProjectPositiveConservative( ...
            z, problem.mass, problem.grid.h, solver.projection_tol);
    otherwise
        error('src:solvers:ApplyProx:UnknownProjection', ...
            'Unknown projection backend "%s".', solver.projection_name);
end
info.backend = lower(char(solver.projection_name));
info.mass_error = abs(src.constraints.Mass(rho, problem.grid.h) - problem.mass);
info.min_rho = min(rho);
info.lambda_iterations = 0;
info.max_inner_iterations = 0;
end
