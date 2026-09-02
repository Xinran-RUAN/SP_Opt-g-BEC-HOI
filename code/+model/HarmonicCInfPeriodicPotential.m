function [V, info] = HarmonicCInfPeriodicPotential(x, L, R0, R1)
%HARMONICCINFPERIODICPOTENTIAL Boundary-only C-infinity continuation.
%
% V=x^2/2 on |x|<=R0, transitions with the standard flat C-infinity
% step on R0<|x|<R1, and equals L^2/2 on |x|>=R1.  Consequently the
% periodic extension is constant near +/-L while the harmonic core is
% unchanged exactly.

if ~isreal(x) || any(~isfinite(x(:))) ...
        || ~isscalar(L) || ~isfinite(L) || L <= 0 ...
        || ~isscalar(R0) || ~isfinite(R0) ...
        || ~isscalar(R1) || ~isfinite(R1) ...
        || R0 <= 0 || R0 >= R1 || R1 >= L
    error('model:HarmonicCInfPeriodicPotential:InvalidInput', ...
        'Require finite real x and 0 < R0 < R1 < L.');
end

harmonic = 0.5 * x .^ 2;
flatValue = 0.5 * L ^ 2;
t = (abs(x) - R0) / (R1 - R0);
chi = model.SmoothStepCInf(t);
V = (1 - chi) .* harmonic + chi .* flatValue;

if nargout > 1
    info.L = L;
    info.R0 = R0;
    info.R1 = R1;
    info.flat_value = flatValue;
    info.chi = chi;
    info.core_mask = abs(x) <= R0;
    info.transition_mask = abs(x) > R0 & abs(x) < R1;
    info.flat_mask = abs(x) >= R1;
    if any(info.core_mask(:))
        info.max_core_error = max(abs( ...
            V(info.core_mask) - harmonic(info.core_mask)));
    else
        info.max_core_error = NaN;
    end
    info.minimum = min(V(:));
end
end
