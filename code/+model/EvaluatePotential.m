function [values, potential] = EvaluatePotential(grid, potential)
%EVALUATEPOTENTIAL Evaluate an inline trapping-potential expression.
%
% Preferred interface:
%   potential.V     = @(x) 0.5*x.^2;
%   potential.label = 'V(x) = x^2/2';

if nargin < 2 || ~isstruct(potential) ...
        || ~isfield(potential, 'V') ...
        || ~isa(potential.V, 'function_handle')
    error('model:EvaluatePotential:InvalidExpression', ...
        'trapping_potential.V must be a function handle of x.');
end
if ~isstruct(grid) || ~all(isfield(grid, {'x', 'N'}))
    error('model:EvaluatePotential:InvalidGrid', ...
        'grid must contain x and N.');
end

values = potential.V(grid.x);
if isscalar(values)
    values = repmat(values, grid.N, 1);
else
    values = values(:);
end
if numel(values) ~= grid.N || ~isreal(values) ...
        || any(~isfinite(values))
    error('model:EvaluatePotential:InvalidValues', ...
        'V(x) must return a finite real scalar or one value per grid node.');
end

if ~isfield(potential, 'label') || isempty(potential.label)
    potential.label = func2str(potential.V);
end
potential.label = char(potential.label);
potential.values = values;
potential.minimum = min(values);
potential.maximum = max(values);
potential.source = 'inline_function_handle';
end
