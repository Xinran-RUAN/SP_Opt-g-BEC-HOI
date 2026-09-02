function fisher = NormalizeFisher(fisher)
%NORMALIZEFISHER Canonicalize the Fisher denominator handle interface.
%
% The production interface is r_epsilon, dr_epsilon, d2r_epsilon.  The
% former s_epsilon spelling is accepted only as a compatibility input for
% archived runs; aliases are attached so old result readers keep working.

canonical = {'r_epsilon', 'dr_epsilon', 'd2r_epsilon'};
legacy = {'s_epsilon', 'ds_epsilon', 'd2s_epsilon'};
if ~isstruct(fisher)
    error('src:regularization:NormalizeFisher:InvalidInput', ...
        'Fisher regularization must be a struct of function handles.');
end
hasCanonical = all(isfield(fisher, canonical));
hasLegacy = all(isfield(fisher, legacy));
if ~hasCanonical && hasLegacy
    fisher.r_epsilon = fisher.s_epsilon;
    fisher.dr_epsilon = fisher.ds_epsilon;
    fisher.d2r_epsilon = fisher.d2s_epsilon;
elseif ~hasCanonical
    error('src:regularization:NormalizeFisher:MissingHandles', ...
        ['Provide r_epsilon, dr_epsilon, and d2r_epsilon ' ...
        'function handles.']);
end
if ~all(cellfun(@(name) isa(fisher.(name), 'function_handle'), canonical))
    error('src:regularization:NormalizeFisher:InvalidHandles', ...
        'The canonical Fisher fields must be function handles.');
end

% Read-only backward aliases.  Production Energy/Gradient/Hessian code
% uses only the canonical fields above.
fisher.s_epsilon = fisher.r_epsilon;
fisher.ds_epsilon = fisher.dr_epsilon;
fisher.d2s_epsilon = fisher.d2r_epsilon;
end
