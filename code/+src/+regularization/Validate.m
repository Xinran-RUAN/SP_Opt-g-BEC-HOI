function [reg, info] = Validate(reg)
%VALIDATE Normalize and validate an active denominator regularization.

if ~isstruct(reg) || ~isfield(reg, 'name') || ~isfield(reg, 'epsilon')
    error('src:regularization:Validate:InvalidInput', ...
        'reg must contain name and epsilon.');
end
reg.name = lower(char(reg.name));
if ~ismember(reg.name, src.regularization.SupportedNames())
    error('src:regularization:Validate:UnsupportedName', ...
        ['Regularization "%s" is not active. Only convexity-preserving ' ...
         'regularizations are accepted.'], reg.name);
end
if ~isscalar(reg.epsilon) || ~isfinite(reg.epsilon) || reg.epsilon <= 0
    error('src:regularization:Validate:InvalidEpsilon', ...
        'epsilon must be a positive finite scalar.');
end
if ~isfield(reg, 'transition_width') || isempty(reg.transition_width)
    reg.transition_width = reg.epsilon;
end
if ~isscalar(reg.transition_width) || ~isfinite(reg.transition_width) || ...
        reg.transition_width <= 0
    error('src:regularization:Validate:InvalidTransitionWidth', ...
        'transition_width must be a positive finite scalar.');
end

info.name = reg.name;
info.is_convex_preserving = true;
info.epsilon = reg.epsilon;
info.transition_width = reg.transition_width;
switch reg.name
    case 'shift_smooth'
        info.regularity = 'C-infinity';
    case 'piecewise_c1'
        info.regularity = 'C1';
    case 'piecewise_c2'
        info.regularity = 'C2';
    case 'piecewise_c3'
        info.regularity = 'C3';
end

kineticExpression = src.regularization.KineticExpression(reg);
reg.kinetic_expression = kineticExpression.text;
info.kinetic_expression = kineticExpression.text;

if ~info.is_convex_preserving
    error('src:regularization:Validate:NotConvexPreserving', ...
        'The active solver requires a convexity-preserving regularization.');
end
end
