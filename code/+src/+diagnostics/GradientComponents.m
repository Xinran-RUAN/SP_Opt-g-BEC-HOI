function components = GradientComponents(rho, problem)
%GRADIENTCOMPONENTS Read-only decomposition of the production gradient.
%
% The derivative and adjoint derivative are evaluated through the same
% Fourier pseudospectral helpers used by src.discretization.ps.Gradient.
% Consequently this diagnostic does not assume a separate D^T convention.

rho = rho(:);
if numel(rho) ~= problem.plan.N
    error('src:diagnostics:GradientComponents:SizeMismatch', ...
        'rho must contain problem.plan.N nodal values.');
end

q = src.discretization.ps.FirstDerivative(rho, problem.plan);
fisher = src.regularization.ResolveFisher(problem);
[r, dr] = src.regularization.EvaluateFisher(rho, fisher);
if any(r <= 0)
    error('src:diagnostics:GradientComponents:NonpositiveFisher', ...
        'r_epsilon(rho) must be strictly positive.');
end

dtQOverR = src.discretization.ps.AdjointFirstDerivative( ...
    q ./ r, problem.plan);
dtQ = src.discretization.ps.AdjointFirstDerivative(q, problem.plan);
[~, dp] = src.potential.Evaluate( ...
    rho, potentialRegularization(problem), fisher.epsilon);

components.fisher = 0.25 * dtQOverR ...
    - 0.125 * (q .^ 2) .* dr ./ (r .^ 2);
components.potential = problem.V(:) .* dp;
components.beta = problem.beta * rho;
components.delta = problem.delta * dtQ;
components.total = components.fisher + components.potential ...
    + components.beta + components.delta;
end

function potential = potentialRegularization(problem)
if isfield(problem, 'potential_regularization')
    potential = problem.potential_regularization;
else
    potential = src.potential.MakeLinear();
end
end
