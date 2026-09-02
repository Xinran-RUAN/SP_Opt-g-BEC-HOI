function diagnostic = EntropyLayerResolution(rho, grid, options)
%ENTROPYLAYERRESOLUTION Conservative nodal resolution of logarithmic tails.
%
% Counts nodes on each side of the main peak whose density lies between
% two configured fractions of the peak. This is a resolution diagnostic,
% not a modification or fitted model of the density.

if nargin < 3
    options = struct();
end
if ~isfield(options, 'lower_relative_density') ...
        || isempty(options.lower_relative_density)
    options.lower_relative_density = 1e-10;
end
if ~isfield(options, 'upper_relative_density') ...
        || isempty(options.upper_relative_density)
    options.upper_relative_density = 1e-2;
end
lowerFraction = options.lower_relative_density;
upperFraction = options.upper_relative_density;
if lowerFraction <= 0 || upperFraction <= lowerFraction ...
        || upperFraction >= 1
    error('src:diagnostics:EntropyLayerResolution:InvalidThresholds', ...
        'Require 0<lower_relative_density<upper_relative_density<1.');
end
rho = rho(:);
if numel(rho) ~= grid.N || any(rho < 0) || any(~isfinite(rho))
    error('src:diagnostics:EntropyLayerResolution:InvalidDensity', ...
        'rho must be a finite nonnegative grid vector.');
end
[peakDensity, peakIndex] = max(rho);
lowerThreshold = lowerFraction * peakDensity;
upperThreshold = upperFraction * peakDensity;
inLayer = rho >= lowerThreshold & rho <= upperThreshold;
leftCount = nnz(inLayer(1:peakIndex));
rightCount = nnz(inLayer(peakIndex:end));

diagnostic.peak_density = peakDensity;
diagnostic.peak_index = peakIndex;
diagnostic.lower_relative_density = lowerFraction;
diagnostic.upper_relative_density = upperFraction;
diagnostic.lower_density_threshold = lowerThreshold;
diagnostic.upper_density_threshold = upperThreshold;
diagnostic.transition_cells_left = leftCount;
diagnostic.transition_cells_right = rightCount;
diagnostic.transition_cells = min(leftCount, rightCount);
diagnostic.transition_width_left = leftCount * grid.h;
diagnostic.transition_width_right = rightCount * grid.h;
diagnostic.layer_node_count = nnz(inLayer);
end
