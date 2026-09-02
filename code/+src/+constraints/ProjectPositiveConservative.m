function rho = ProjectPositiveConservative(z, mass, h, tolerance)
%PROJECTPOSITIVECONSERVATIVE Semismooth projection onto the feasible set.
%
% This independently solves h*sum(max(z+lambda,0))=mass by a safeguarded
% semismooth Newton iteration. It does not call ProjectSimplex.

if nargin < 4 || isempty(tolerance)
    tolerance = 1e-13;
end
if ~isscalar(mass) || ~isfinite(mass) || mass <= 0 || ...
        ~isscalar(h) || ~isfinite(h) || h <= 0 || ...
        ~isscalar(tolerance) || tolerance <= 0
    error('src:constraints:ProjectPositiveConservative:InvalidParameters', ...
        'mass, h, and tolerance must be positive finite scalars.');
end
wasRow = isrow(z);
z = z(:);
if isempty(z) || any(~isfinite(z)) || ~isreal(z)
    error('src:constraints:ProjectPositiveConservative:InvalidInput', ...
        'z must be a nonempty finite real vector.');
end

lower = -max(z) - mass / h;
upper = mass / h - min(z);
lambda = min(max(0, lower), upper);

for iteration = 1:100
    shifted = z + lambda;
    positive = shifted > 0;
    residual = h * sum(shifted(positive)) - mass;
    if abs(residual) <= tolerance
        break;
    end
    if residual > 0
        upper = lambda;
    else
        lower = lambda;
    end

    derivative = h * nnz(positive);
    if derivative > 0
        trial = lambda - residual / derivative;
    else
        trial = 0.5 * (lower + upper);
    end
    if ~isfinite(trial) || trial <= lower || trial >= upper
        trial = 0.5 * (lower + upper);
    end
    lambda = trial;
end

rho = max(z + lambda, 0);
% On large tensor grids, forming z+lambda can lose a few ulps when the
% two terms are large and nearly cancel.  A uniform correction on the
% positive set is exactly a final multiplier update for the same
% projection problem; it improves mass summation without clipping or
% changing the active set.
for correctionIteration = 1:3
    positive = rho > 0;
    massResidual = h * sum(rho) - mass;
    if massResidual == 0 || ~any(positive)
        break;
    end
    correction = massResidual / (h * nnz(positive));
    corrected = rho;
    corrected(positive) = corrected(positive) - correction;
    if any(corrected(positive) <= 0)
        break;
    end
    rho = corrected;
end
massError = abs(h * sum(rho) - mass);
if massError > max(10 * tolerance, 500 * eps(max(1, mass)))
    error('src:constraints:ProjectPositiveConservative:NoConvergence', ...
        'Semismooth projection mass residual is %.3e.', massError);
end
if wasRow
    rho = rho.';
end
end
