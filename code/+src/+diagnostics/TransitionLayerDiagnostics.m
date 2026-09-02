function diagnostic = TransitionLayerDiagnostics( ...
    rho, grid, plan, potential, options)
%TRANSITIONLAYERDIAGNOSTICS Resolve q_sigma=dp_sigma(rho) transitions.
%
% The primary layer is defined only by q_sigma in [0.1,0.9]. Crossing
% locations are linearly interpolated between adjacent collocation nodes.
% Missing crossings remain NaN and carry an explicit reason; they are never
% replaced by a zero width or a density-tolerance surrogate.

if nargin < 5
    options = struct();
end
options = defaults(options);
rho = rho(:);
if numel(rho) ~= grid.N || plan.N ~= grid.N
    error('src:diagnostics:TransitionLayerDiagnostics:SizeMismatch', ...
        'rho, grid, and plan sizes must agree.');
end
required = {'dp_sigma', 'd2p_sigma'};
if ~isstruct(potential) || ~all(isfield(potential, required))
    error('src:diagnostics:TransitionLayerDiagnostics:MissingHandles', ...
        'potential must provide dp_sigma and d2p_sigma.');
end

qSigma = potential.dp_sigma(rho);
d2p = potential.d2p_sigma(rho);
if any(~isfinite(qSigma)) || any(~isfinite(d2p))
    error('src:diagnostics:TransitionLayerDiagnostics:NonfiniteValues', ...
        'Potential derivative handles returned nonfinite values.');
end

drho = src.discretization.ps.FirstDerivative(rho, plan);
d2rho = src.discretization.ps.FirstDerivative(drho, plan);
d3rho = src.discretization.ps.FirstDerivative(d2rho, plan);
dqdx = d2p .* drho;

primary = measureBand(grid.x, qSigma, grid.L, grid.h, ...
    options.primary_levels(1), options.primary_levels(2));
secondary = measureBand(grid.x, qSigma, grid.L, grid.h, ...
    options.secondary_levels(1), options.secondary_levels(2));
primary = localScales(primary, qSigma, dqdx, drho, d2rho, d3rho, d2p);
secondary = localScales(secondary, qSigma, dqdx, drho, d2rho, d3rho, d2p);

diagnostic.sigma = getFieldOr(potential, 'sigma', NaN);
diagnostic.q_sigma = qSigma;
diagnostic.dqdx = dqdx;
diagnostic.drho = drho;
diagnostic.d2rho = d2rho;
diagnostic.d3rho = d3rho;
diagnostic.primary = primary;
diagnostic.secondary = secondary;
diagnostic.width_mean = primary.width_mean;
diagnostic.points_per_layer = primary.points_per_layer;
diagnostic.n_nodes_layer = primary.n_nodes_layer;
diagnostic.x_layer_right = primary.center_right;
diagnostic.x_layer_left = primary.center_left;
diagnostic.max_abs_dqdx = max(abs(dqdx));
if diagnostic.max_abs_dqdx > 0
    diagnostic.width_slope = diff(options.primary_levels) ...
        / diagnostic.max_abs_dqdx;
    diagnostic.points_per_slope_width = ...
        diagnostic.width_slope / grid.h;
else
    diagnostic.width_slope = Inf;
    diagnostic.points_per_slope_width = Inf;
end
diagnostic.resolution_label = resolutionLabel(primary.points_per_layer);
diagnostic.slope_resolution_label = ...
    resolutionLabel(diagnostic.points_per_slope_width);
diagnostic.min_density = min(rho);
diagnostic.min_q_sigma = min(qSigma);
diagnostic.max_q_sigma = max(qSigma);
end

function band = measureBand(x, q, L, h, lowLevel, highLevel)
[xHighRight, highRightReason] = crossingOnSide( ...
    x, q, highLevel, 'right');
[xLowRight, lowRightReason] = crossingOnSide( ...
    x, q, lowLevel, 'right');
[xHighLeft, highLeftReason] = crossingOnSide( ...
    x, q, highLevel, 'left');
[xLowLeft, lowLeftReason] = crossingOnSide( ...
    x, q, lowLevel, 'left');
midLevel = 0.5 * (lowLevel + highLevel);
[xMidRight, midRightReason] = crossingOnSide( ...
    x, q, midLevel, 'right');
[xMidLeft, midLeftReason] = crossingOnSide( ...
    x, q, midLevel, 'left');

band.low_level = lowLevel;
band.high_level = highLevel;
band.x_high_right = xHighRight;
band.x_low_right = xLowRight;
band.x_high_left = xHighLeft;
band.x_low_left = xLowLeft;
band.mid_level = midLevel;
band.x_mid_right = xMidRight;
band.x_mid_left = xMidLeft;
band.crossing_reasons.high_right = highRightReason;
band.crossing_reasons.low_right = lowRightReason;
band.crossing_reasons.high_left = highLeftReason;
band.crossing_reasons.low_left = lowLeftReason;
band.crossing_reasons.mid_right = midRightReason;
band.crossing_reasons.mid_left = midLeftReason;
band.valid = all(isfinite([xHighRight, xLowRight, xHighLeft, xLowLeft]));
band.mask = q >= lowLevel & q <= highLevel;
band.n_nodes_layer = nnz(band.mask);

if band.valid
    band.width_right = abs(xLowRight - xHighRight);
    band.width_left = abs(xHighLeft - xLowLeft);
    band.width_mean = 0.5 * (band.width_left + band.width_right);
    band.points_per_layer = band.width_mean / h;
    band.asymmetry = abs(band.width_left - band.width_right) ...
        / max(band.width_mean, eps);
    band.center_right = 0.5 * (xLowRight + xHighRight);
    band.center_left = 0.5 * (xLowLeft + xHighLeft);
    band.reason = 'all crossings found';
else
    band.width_right = NaN;
    band.width_left = NaN;
    band.width_mean = NaN;
    band.points_per_layer = NaN;
    band.asymmetry = NaN;
    band.center_right = NaN;
    band.center_left = NaN;
    band.reason = joinReasons(band.crossing_reasons);
end

if isfinite(xHighRight) && isfinite(xHighLeft)
    rightLowerBound = L - xHighRight;
    leftLowerBound = xHighLeft + L;
    band.censored_width_lower_bound = ...
        0.5 * (rightLowerBound + leftLowerBound);
    band.censored_points_lower_bound = ...
        band.censored_width_lower_bound / h;
else
    band.censored_width_lower_bound = NaN;
    band.censored_points_lower_bound = NaN;
end
end

function band = localScales(band, q, dqdx, drho, d2rho, d3rho, d2p)
mask = band.mask;
if any(mask)
    band.max_abs_dqdx = max(abs(dqdx(mask)));
    band.max_abs_drho = max(abs(drho(mask)));
    band.max_abs_d2rho = max(abs(d2rho(mask)));
    band.max_abs_d3rho = max(abs(d3rho(mask)));
    band.max_d2p = max(d2p(mask));
    band.max_abs_d2p_drho = max(abs(d2p(mask) .* drho(mask)));
    band.q_min_on_nodes = min(q(mask));
    band.q_max_on_nodes = max(q(mask));
else
    names = {'max_abs_dqdx', 'max_abs_drho', 'max_abs_d2rho', ...
        'max_abs_d3rho', 'max_d2p', 'max_abs_d2p_drho', ...
        'q_min_on_nodes', 'q_max_on_nodes'};
    for j = 1:numel(names)
        band.(names{j}) = NaN;
    end
end
end

function [location, reason] = crossingOnSide(x, q, level, side)
switch lower(side)
    case 'right'
        indices = find(x >= 0);
    case 'left'
        indices = flipud(find(x <= 0));
    otherwise
        error('Unknown side "%s".', side);
end
xSide = x(indices);
qSide = q(indices);
crossing = find(qSide(1:end-1) >= level ...
    & qSide(2:end) <= level, 1, 'first');
if isempty(crossing)
    location = NaN;
    if max(qSide) < level
        reason = sprintf('q never reaches %.3g on %s side', level, side);
    elseif min(qSide) > level
        reason = sprintf('q never falls to %.3g on %s side', level, side);
    else
        reason = sprintf('no outward %.3g crossing on %s side', level, side);
    end
    return;
end
x1 = xSide(crossing);
x2 = xSide(crossing + 1);
q1 = qSide(crossing);
q2 = qSide(crossing + 1);
if q2 == q1
    location = 0.5 * (x1 + x2);
else
    location = x1 + (level - q1) * (x2 - x1) / (q2 - q1);
end
reason = 'crossing found';
end

function text = joinReasons(reasons)
names = fieldnames(reasons);
parts = {};
for j = 1:numel(names)
    if ~strcmp(reasons.(names{j}), 'crossing found')
        parts{end + 1} = reasons.(names{j}); %#ok<AGROW>
    end
end
if isempty(parts)
    text = 'crossing data incomplete';
else
    text = strjoin(parts, '; ');
end
end

function label = resolutionLabel(points)
if ~isfinite(points)
    label = 'not measurable';
elseif points < 2
    label = 'severely unresolved';
elseif points < 6
    label = 'marginally resolved';
elseif points < 12
    label = 'moderately resolved';
else
    label = 'well resolved';
end
end

function options = defaults(options)
if ~isfield(options, 'primary_levels') || isempty(options.primary_levels)
    options.primary_levels = [0.1, 0.9];
end
if ~isfield(options, 'secondary_levels') || isempty(options.secondary_levels)
    options.secondary_levels = [0.25, 0.75];
end
if any(diff(options.primary_levels) <= 0) ...
        || any(diff(options.secondary_levels) <= 0)
    error('Transition levels must be strictly increasing [low, high].');
end
end

function value = getFieldOr(data, name, fallback)
if isfield(data, name)
    value = data.(name);
else
    value = fallback;
end
end
