function fisher = ResolveFisher(problem)
%RESOLVEFISHER Return the canonical r_epsilon handle interface.

if isfield(problem, 'fisher_regularization') ...
        && ~isempty(problem.fisher_regularization)
    fisher = problem.fisher_regularization;
elseif isfield(problem, 'regularization')
    fisher = src.regularization.MakeBuiltIn(problem.regularization);
else
    error('src:regularization:ResolveFisher:MissingInterface', ...
        'problem must contain fisher_regularization or legacy regularization.');
end
fisher = src.regularization.NormalizeFisher(fisher);
end
