function diagnostic = FourierEnvelopeFit(values, options)
%FOURIERENVELOPEFIT Band envelopes and two descriptive decay fits.
%
% Fits are diagnostics only. Coefficients at the floating-point noise floor
% are excluded so a finite-band trigonometric polynomial is not classified
% from an arbitrary regression through FFT roundoff.

if nargin < 2
    options = struct();
end
options = defaults(options);
[coefficients, modes] = src.discretization.ps.FourierCoefficients(values);
amplitudes = abs(coefficients);
absoluteModes = abs(modes);

bands = options.bands;
bandCount = size(bands, 1);
bandMax = NaN(bandCount, 1);
bandRms = NaN(bandCount, 1);
bandSamples = zeros(bandCount, 1);
for j = 1:bandCount
    mask = absoluteModes >= bands(j, 1) ...
        & absoluteModes <= bands(j, 2);
    bandSamples(j) = nnz(mask);
    if any(mask)
        bandMax(j) = max(amplitudes(mask));
        bandRms(j) = sqrt(mean(amplitudes(mask) .^ 2));
    end
end

scale = max(amplitudes);
noiseFloor = max(options.absolute_floor, options.relative_floor * scale);
fitMask = absoluteModes >= options.fit_range(1) ...
    & absoluteModes <= min(options.fit_range(2), numel(values) / 2 - 1) ...
    & amplitudes > noiseFloor;
k = absoluteModes(fitMask);
logAmplitude = log(amplitudes(fitMask));

[algebraicSlope, r2Algebraic] = fitLine(log(k), logAmplitude);
[exponentialSlope, r2Exponential] = fitLine(k, logAmplitude);

diagnostic.modes = modes;
diagnostic.coefficients = coefficients;
diagnostic.amplitudes = amplitudes;
diagnostic.bands = bands;
diagnostic.band_max = bandMax;
diagnostic.band_rms = bandRms;
diagnostic.band_sample_count = bandSamples;
diagnostic.fit_range = options.fit_range;
diagnostic.fit_count = nnz(fitMask);
diagnostic.noise_floor = noiseFloor;
diagnostic.algebraic_slope = -algebraicSlope;
diagnostic.R2_algebraic = r2Algebraic;
diagnostic.exponential_slope = -exponentialSlope;
diagnostic.R2_exponential = r2Exponential;
diagnostic.fit_modes = k;
end

function options = defaults(options)
if ~isfield(options, 'bands') || isempty(options.bands)
    options.bands = [8, 16; 16, 32; 32, 64; 64, 128; ...
        128, 256; 256, 512];
end
if ~isfield(options, 'fit_range') || isempty(options.fit_range)
    options.fit_range = [8, 512];
end
if ~isfield(options, 'relative_floor') || isempty(options.relative_floor)
    options.relative_floor = 100 * eps;
end
if ~isfield(options, 'absolute_floor') || isempty(options.absolute_floor)
    options.absolute_floor = realmin;
end
if ~isequal(size(options.fit_range), [1, 2]) ...
        || options.fit_range(1) <= 0 ...
        || options.fit_range(2) < options.fit_range(1)
    error('src:diagnostics:FourierEnvelopeFit:InvalidFitRange', ...
        'fit_range must be [positiveMinimum, maximum].');
end
if size(options.bands, 2) ~= 2 || any(options.bands(:) < 0)
    error('src:diagnostics:FourierEnvelopeFit:InvalidBands', ...
        'bands must contain nonnegative [lower, upper] rows.');
end
end

function [slope, rSquared] = fitLine(abscissa, ordinate)
if numel(abscissa) < 4 || max(abscissa) == min(abscissa)
    slope = NaN;
    rSquared = NaN;
    return;
end
coefficients = polyfit(abscissa, ordinate, 1);
prediction = polyval(coefficients, abscissa);
residualPower = sum((ordinate - prediction) .^ 2);
totalPower = sum((ordinate - mean(ordinate)) .^ 2);
slope = coefficients(1);
if totalPower == 0
    rSquared = NaN;
else
    rSquared = 1 - residualPower / totalPower;
end
end
