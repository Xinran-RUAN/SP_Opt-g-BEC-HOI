function [fisher, info] = ValidateFisher(fisher, rhoTest)
%VALIDATEFISHER Validate the canonical r_epsilon handle interface.

if nargin < 2 || isempty(rhoTest)
    rhoTest = [0; 1];
end
fisher = src.regularization.NormalizeFisher(fisher);
if ~isfield(fisher, 'epsilon') || ~isscalar(fisher.epsilon) ...
        || ~isfinite(fisher.epsilon) || fisher.epsilon <= 0
    error('src:regularization:ValidateFisher:InvalidEpsilon', ...
        'fisher_regularization.epsilon must be positive and finite.');
end
if ~isfield(fisher, 'label') || isempty(fisher.label)
    fisher.label = func2str(fisher.r_epsilon);
end
[r, dr, d2r] = src.regularization.EvaluateFisher(rhoTest, fisher);
scale = max([1; abs(r(:)); abs(dr(:)); abs(d2r(:))]);
tolerance = 1e3 * eps(scale);
if any(r(:) <= 0)
    error('src:regularization:ValidateFisher:NonpositiveR', ...
        'r_epsilon must be strictly positive on the tested nodal range.');
end
if any(dr(:) < -tolerance)
    error('src:regularization:ValidateFisher:NegativeDerivative', ...
        'dr_epsilon must be nonnegative for the convex solver.');
end
if any(d2r(:) > tolerance)
    error('src:regularization:ValidateFisher:PositiveCurvature', ...
        'd2r_epsilon must be nonpositive for convexity preservation.');
end
fisher.is_convex_preserving = true;
info.epsilon = fisher.epsilon;
info.label = fisher.label;
info.r_epsilon = fisher.r_epsilon;
info.dr_epsilon = fisher.dr_epsilon;
info.d2r_epsilon = fisher.d2r_epsilon;
% Compatibility aliases retained in MAT archives created by older runs.
info.s_epsilon = fisher.r_epsilon;
info.ds_epsilon = fisher.dr_epsilon;
info.d2s_epsilon = fisher.d2r_epsilon;
info.is_convex_preserving = true;
info.minimum_r_test = min(r(:));
info.minimum_dr_test = min(dr(:));
info.maximum_d2r_test = max(d2r(:));
end
