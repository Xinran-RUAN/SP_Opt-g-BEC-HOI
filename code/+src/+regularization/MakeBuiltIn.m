function fisher = MakeBuiltIn(reg)
%MAKEBUILTIN Adapt a legacy named Fisher regularization to paper handles.

[reg, info] = src.regularization.Validate(reg);
epsilon = reg.epsilon;
transitionWidth = reg.transition_width;
name = reg.name;
switch name
    case 'shift_smooth'
        fisher.r_epsilon = @(rho) rho + epsilon;
        fisher.dr_epsilon = @(rho) ones(size(rho));
        fisher.d2r_epsilon = @(rho) zeros(size(rho));
        label = 'r_epsilon(rho) = rho + epsilon';
    otherwise
        continuity = sscanf(name, 'piecewise_c%d');
        if ismember(continuity, [1, 2])
            % Formal C1/C2 family: the matching threshold equals epsilon.
            fisher = src.regularization.PiecewiseFisherCm( ...
                epsilon, continuity);
            label = fisher.label;
        else
            % Retain the archived C3 construction only in this legacy
            % adapter; it is not part of the current formal experiment.
            fisher.r_epsilon = @(rho) piecewiseField( ...
                rho, epsilon, transitionWidth, continuity, 0);
            fisher.dr_epsilon = @(rho) piecewiseField( ...
                rho, epsilon, transitionWidth, continuity, 1);
            fisher.d2r_epsilon = @(rho) piecewiseField( ...
                rho, epsilon, transitionWidth, continuity, 2);
            label = sprintf('legacy r_epsilon: %s, width %.6g', ...
                name, transitionWidth);
        end
end
fisher.epsilon = epsilon;
fisher.label = label;
fisher.source = 'legacy_builtin';
fisher.legacy_name = name;
fisher.is_convex_preserving = info.is_convex_preserving;
fisher = src.regularization.NormalizeFisher(fisher);
end

function value = piecewiseField(rho, epsilon, width, continuity, order)
roundoff = 100 * eps(max(1, max(abs(rho(:)))));
if any(rho(:) < -roundoff)
    error('src:regularization:MakeBuiltIn:NegativeDensity', ...
        'The legacy piecewise Fisher map is defined for rho >= 0.');
end
r = max(rho, 0);
inside = r <= width;
z = 1 - r(inside) / width;
switch order
    case 0
        value = epsilon + r + width / (continuity + 1);
        value(inside) = epsilon + r(inside) + width / (continuity + 1) ...
            .* (1 - z .^ (continuity + 1));
    case 1
        value = ones(size(r));
        value(inside) = 1 + z .^ continuity;
    case 2
        value = zeros(size(r));
        value(inside) = -(continuity / width) ...
            .* z .^ (continuity - 1);
end
end
