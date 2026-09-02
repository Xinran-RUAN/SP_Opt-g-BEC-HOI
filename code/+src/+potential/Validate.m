function [potential, info] = Validate(potential, epsilon, rhoTest)
%VALIDATE Normalize p_sigma handles and validate convexity assumptions.
%
% Inline configurations already containing p_sigma, dp_sigma, and
% d2p_sigma are used directly. Legacy names are adapted once here.

if nargin < 1 || isempty(potential)
    potential = src.potential.MakeLinear();
end
if nargin < 2 || ~isscalar(epsilon) || ~isfinite(epsilon) || epsilon <= 0
    error('src:potential:Validate:InvalidEpsilon', ...
        'epsilon must be a positive finite scalar.');
end
if nargin < 3 || isempty(rhoTest)
    rhoTest = [0; epsilon; 1];
end
required = {'p_sigma', 'dp_sigma', 'd2p_sigma'};
hasHandles = isstruct(potential) && all(isfield(potential, required)) ...
    && all(cellfun(@(field) isa(potential.(field), 'function_handle'), required));
if ~hasHandles
    potential = adaptLegacy(potential, epsilon);
end
if ~isfield(potential, 'sigma') || ~isscalar(potential.sigma) ...
        || ~isfinite(potential.sigma) || potential.sigma < 0
    error('src:potential:Validate:InvalidSigma', ...
        'potential_regularization.sigma must be finite and nonnegative.');
end
if ~isfield(potential, 'label') || isempty(potential.label)
    potential.label = func2str(potential.p_sigma);
end
if ~isfield(potential, 'name') || isempty(potential.name)
    potential.name = 'custom_inline';
end
if ~isfield(potential, 'power') || isempty(potential.power)
    potential.power = NaN;
end
if ~isfield(potential, 'prox_type') || isempty(potential.prox_type)
    potential.prox_type = 'generic_convex';
end

rhoTest = max(rhoTest(:), 0);
p = potential.p_sigma(rhoTest);
dp = potential.dp_sigma(rhoTest);
d2p = potential.d2p_sigma(rhoTest);
validateOutputs(p, dp, d2p, rhoTest);
valueTolerance = 1e3 * eps(max([1; abs(p(:))]));
slopeTolerance = 1e3 * eps(max([1; abs(dp(:))]));
curvatureTolerance = 1e3 * eps(max([1; abs(d2p(:))]));
p0 = potential.p_sigma(0);
if ~isscalar(p0) || ~isfinite(p0) || abs(p0) > valueTolerance
    error('src:potential:Validate:NonzeroVacuumValue', ...
        'p_sigma(0) must equal zero to roundoff.');
end
if any(dp(:) < -slopeTolerance)
    error('src:potential:Validate:NegativeSlope', ...
        'dp_sigma must be nonnegative on the tested nodal range.');
end
if any(d2p(:) < -curvatureTolerance)
    error('src:potential:Validate:NonconvexPotential', ...
        'd2p_sigma must be nonnegative on the tested nodal range.');
end

info.name = char(potential.name);
info.family = ternary(strcmpi(potential.prox_type, 'linear'), ...
    'linear', 'generic_convex');
info.power = potential.power;
info.sigma = potential.sigma;
info.label = potential.label;
info.p_sigma = potential.p_sigma;
info.dp_sigma = potential.dp_sigma;
info.d2p_sigma = potential.d2p_sigma;
info.prox_type = potential.prox_type;
info.is_convex = true;
info.dp_at_zero = potential.dp_sigma(0);
info.minimum_dp_test = min(dp(:));
info.minimum_d2p_test = min(d2p(:));
end

function potential = adaptLegacy(potential, epsilon)
if ~isstruct(potential) || ~isfield(potential, 'name') ...
        || isempty(potential.name)
    potential = src.potential.MakeLinear();
    return;
end
name = lower(char(potential.name));
if ~ismember(name, src.potential.SupportedNames())
    error('src:potential:Validate:UnsupportedName', ...
        'Unsupported potential regularization "%s".', name);
end
if strcmp(name, 'linear')
    potential = src.potential.MakeLinear();
    return;
end
switch name
    case 'sqrt_same_scale'
        power = 1;
    case 'sqrt_squared_scale'
        power = 2;
    case 'sqrt_power'
        if ~isfield(potential, 'power') || isempty(potential.power)
            potential.power = 2;
        end
        power = potential.power;
end
if ~isscalar(power) || ~isfinite(power) || power < 1 || power > 4 ...
        || power ~= round(power)
    error('src:potential:Validate:InvalidPower', ...
        'sqrt_power requires an integer power in [1,4].');
end
sigma = epsilon ^ power;
potential.sigma = sigma;
potential.p_sigma = @(rho) ...
    rho .^ 2 ./ (hypot(rho, sigma) + sigma);
potential.dp_sigma = @(rho) rho ./ hypot(rho, sigma);
potential.d2p_sigma = @(rho) ...
    (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
potential.label = 'p_sigma(rho) = sqrt(rho^2 + sigma^2) - sigma';
potential.prox_type = 'generic_convex';
potential.name = name;
potential.power = power;
end

function validateOutputs(p, dp, d2p, rho)
if ~isequal(size(p), size(rho)) || ~isequal(size(dp), size(rho)) ...
        || ~isequal(size(d2p), size(rho))
    error('src:potential:Validate:SizeMismatch', ...
        'p_sigma and its derivatives must preserve the size of rho.');
end
if ~isreal(p) || ~isreal(dp) || ~isreal(d2p) ...
        || any(~isfinite(p(:))) || any(~isfinite(dp(:))) ...
        || any(~isfinite(d2p(:)))
    error('src:potential:Validate:NonfiniteOutput', ...
        'p_sigma and its derivatives must be finite and real.');
end
end

function value = ternary(condition, trueValue, falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
