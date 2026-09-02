function [energy, parts] = Energy(rho, problem)
%ENERGY Discrete regularized HOI density energy on a periodic grid.

rho = rho(:);
q = src.discretization.ps.SpatialGradient(rho, problem.plan);
qSquared = sum(q .^ 2, 2);
[r, ~, ~] = src.regularization.EvaluateFisher( ...
    rho, src.regularization.ResolveFisher(problem));
[potentialDensity, ~] = src.potential.Evaluate( ...
    rho, potentialRegularization(problem), fisherEpsilon(problem));

parts.fisher = problem.grid.h * sum(qSquared ./ (8 .* r));
parts.kinetic = parts.fisher; % compatibility alias
parts.potential = problem.grid.h * sum(problem.V(:) .* potentialDensity);
parts.beta = problem.grid.h * (problem.beta / 2) * sum(rho .^ 2);
parts.delta = problem.grid.h * (problem.delta / 2) * sum(qSquared);
energy = parts.kinetic + parts.potential + parts.beta + parts.delta;
end

function epsilon = fisherEpsilon(problem)
fisher = src.regularization.ResolveFisher(problem);
epsilon = fisher.epsilon;
end

function potreg = potentialRegularization(problem)
if isfield(problem, 'potential_regularization')
    potreg = problem.potential_regularization;
else
    potreg = src.potential.MakeLinear();
end
end
