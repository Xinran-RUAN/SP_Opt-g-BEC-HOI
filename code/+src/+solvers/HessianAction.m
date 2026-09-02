function hessianEta = HessianAction(rho, eta, modelData, plan, reg)
%HESSIANACTION Matrix-free action of the discrete variational Hessian.
%
% For piecewise C1, d2r_epsilon is the selected one-sided/generalized
% second derivative supplied by the legacy adapter. It is used by the
% semismooth polish and is not a claim of classical C2 regularity.

rho = rho(:);
eta = eta(:);
if numel(rho) ~= plan.N || numel(eta) ~= plan.N
    error('src:solvers:HessianAction:SizeMismatch', ...
        'rho and eta must have plan.N entries.');
end
q = src.discretization.ps.SpatialGradient(rho, plan);
dq = src.discretization.ps.SpatialGradient(eta, plan);
if isfield(modelData, 'fisher_regularization')
    fisher = modelData.fisher_regularization;
else
    fisher = src.regularization.MakeBuiltIn(reg);
end
[r, dr, d2r] = src.regularization.EvaluateFisher(rho, fisher);

insideAdjoint = dq ./ r - q .* (dr .* eta ./ (r .^ 2));
firstTerm = 0.25 * src.discretization.ps.AdjointSpatialGradient( ...
    insideAdjoint, plan);
denominatorDerivative = d2r ./ (r .^ 2) ...
    - 2 * (dr .^ 2) ./ (r .^ 3);
secondTerm = -0.125 * ( ...
    2 * sum(q .* dq, 2) .* dr ./ (r .^ 2) ...
    + sum(q .^ 2, 2) .* denominatorDerivative .* eta);
dtDeta = src.discretization.ps.AdjointSpatialGradient(dq, plan);
hessianEta = firstTerm + secondTerm ...
    + modelData.beta * eta + modelData.delta * dtDeta;
end
