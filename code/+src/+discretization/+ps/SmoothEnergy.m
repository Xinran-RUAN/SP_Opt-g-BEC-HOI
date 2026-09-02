function [energy, parts] = SmoothEnergy(rho, problem)
%SMOOTHENERGY Fisher, beta, and delta part of the target energy.
%
% The potential contribution is deliberately excluded. r_epsilon is read
% from the normalized Fisher handle interface.

rho = rho(:);
q = src.discretization.ps.SpatialGradient(rho, problem.plan);
qSquared = sum(q .^ 2, 2);
[r, ~, ~] = src.regularization.EvaluateFisher( ...
    rho, src.regularization.ResolveFisher(problem));

parts.fisher = problem.grid.h * sum(qSquared ./ (8 .* r));
parts.kinetic = parts.fisher; % compatibility alias
parts.beta = problem.grid.h * (problem.beta / 2) * sum(rho .^ 2);
parts.delta = problem.grid.h * (problem.delta / 2) * sum(qSquared);
energy = parts.fisher + parts.beta + parts.delta;
end
