function hessianEta = FullHessianAction(rho, eta, problem)
%FULLHESSIANACTION Matrix-free Hessian action of the complete target.

rho = rho(:);
eta = eta(:);
if numel(rho) ~= problem.grid.N || numel(eta) ~= problem.grid.N
    error('src:solvers:FullHessianAction:SizeMismatch', ...
        'rho and eta must have grid.N entries.');
end
hessianEta = src.solvers.HessianAction( ...
    rho, eta, problem, problem.plan, problem.regularization);
potreg = problem.potential_regularization;
fisher = src.regularization.ResolveFisher(problem);
[~, ~, potentialSecondDerivative, ~] = src.potential.Evaluate( ...
    rho, potreg, fisher.epsilon);
potentialCurvature = problem.V(:) .* potentialSecondDerivative;
hessianEta = hessianEta + potentialCurvature .* eta;
end
