function [r, dr, d2r] = EvaluateFisher(rho, fisher)
%EVALUATEFISHER Evaluate r_epsilon and its first two derivatives.

fisher = src.regularization.NormalizeFisher(fisher);
r = fisher.r_epsilon(rho);
dr = fisher.dr_epsilon(rho);
d2r = fisher.d2r_epsilon(rho);
if ~isequal(size(r), size(rho)) || ~isequal(size(dr), size(rho)) ...
        || ~isequal(size(d2r), size(rho))
    error('src:regularization:EvaluateFisher:SizeMismatch', ...
        'Fisher handles must preserve the size of rho.');
end
if ~isreal(r) || ~isreal(dr) || ~isreal(d2r) ...
        || any(~isfinite(r(:))) || any(~isfinite(dr(:))) ...
        || any(~isfinite(d2r(:)))
    error('src:regularization:EvaluateFisher:NonfiniteOutput', ...
        'r_epsilon and its derivatives must be finite and real.');
end
end
