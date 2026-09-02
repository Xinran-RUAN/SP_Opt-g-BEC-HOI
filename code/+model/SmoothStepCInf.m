function chi = SmoothStepCInf(t)
%SMOOTHSTEPCINF Standard flat C-infinity step from zero to one.
%
% chi(t)=0 for t<=0, chi(t)=1 for t>=1, and on 0<t<1
% chi=exp(-1/t)/(exp(-1/t)+exp(-1/(1-t))). Logical masks ensure that
% neither reciprocal is evaluated at an endpoint.

if ~isreal(t) || any(~isfinite(t(:)))
    error('model:SmoothStepCInf:InvalidInput', ...
        't must contain finite real values.');
end
chi = zeros(size(t), 'like', t);
chi(t >= 1) = 1;
interior = t > 0 & t < 1;
if any(interior(:))
    local = t(interior);
    left = exp(-1 ./ local);
    right = exp(-1 ./ (1 - local));
    chi(interior) = left ./ (left + right);
end
end
