function [p, dp, varargout] = Evaluate(rho, potential, epsilon)
%EVALUATE Evaluate p_sigma and its first two derivatives from handles.
%
% Backward-compatible three-output form:
%   [p, dp, info] = Evaluate(...)
% Curvature-aware four-output form:
%   [p, dp, d2p, info] = Evaluate(...)

[potential, info] = src.potential.Validate(potential, epsilon, rho);
p = potential.p_sigma(rho);
dp = potential.dp_sigma(rho);
d2p = potential.d2p_sigma(rho);
if nargout == 3
    info.d2p = d2p;
    varargout{1} = info;
elseif nargout >= 4
    varargout{1} = d2p;
    varargout{2} = info;
end
end
