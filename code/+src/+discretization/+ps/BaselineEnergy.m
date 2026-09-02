function [energy, parts] = BaselineEnergy(rho, problem)
%BASELINEENERGY Evaluate the original linear-potential regularized energy.

baselineProblem = problem;
baselineProblem.potential_regularization = src.potential.MakeLinear();
[energy, parts] = src.discretization.ps.Energy(rho, baselineProblem);
end
