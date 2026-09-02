function gradient = SmoothGradient(rho, problem)
%SMOOTHGRADIENT Variational gradient of Fisher, beta, and delta terms.
%
% The vector follows <u,v>_h=h*sum(u.*v), so no outer factor h appears.
% The potential derivative is deliberately excluded.

rho = rho(:);
q = src.discretization.ps.SpatialGradient(rho, problem.plan);
qSquared = sum(q .^ 2, 2);
[r, dr] = src.regularization.EvaluateFisher( ...
    rho, src.regularization.ResolveFisher(problem));
dtQOverR = src.discretization.ps.AdjointSpatialGradient( ...
    q ./ r, problem.plan);
dtQ = src.discretization.ps.AdjointSpatialGradient(q, problem.plan);

gradient = 0.25 * dtQOverR ...
    - 0.125 * qSquared .* dr ./ (r .^ 2) ...
    + problem.beta * rho + problem.delta * dtQ;
end
