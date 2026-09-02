function gradient = Gradient(rho, problem)
%GRADIENT Variational gradient under <u,v>_h = h*sum(u.*v).
%
% No outer factor h appears in this vector. Energy and gradient share the
% r_epsilon and dr_epsilon are supplied by the normalized Fisher handles.

rho = rho(:);
q = src.discretization.ps.SpatialGradient(rho, problem.plan);
qSquared = sum(q .^ 2, 2);
fisher = src.regularization.ResolveFisher(problem);
[r, dr] = src.regularization.EvaluateFisher(rho, fisher);
dtQOverR = src.discretization.ps.AdjointSpatialGradient(q ./ r, problem.plan);
dtQ = src.discretization.ps.AdjointSpatialGradient(q, problem.plan);
[~, potentialDerivative] = src.potential.Evaluate( ...
    rho, potentialRegularization(problem), fisher.epsilon);

gradient = 0.25 * dtQOverR ...
    - 0.125 * qSquared .* dr ./ (r .^ 2) ...
    + problem.V(:) .* potentialDerivative ...
    + problem.beta * rho + problem.delta * dtQ;
end

function potreg = potentialRegularization(problem)
if isfield(problem, 'potential_regularization')
    potreg = problem.potential_regularization;
else
    potreg = src.potential.MakeLinear();
end
end
