function diagnostic = RelativeSupportDiagnostics(rho, grid)
%RELATIVESUPPORTDIAGNOSTICS Radii of relative-density level crossings.

rho = rho(:);
[peak, peakIndex] = max(rho);
levels = [1e-3, 1e-6, 1e-10, 1e-14];
radii = zeros(size(levels));
leftRadii = zeros(size(levels));
rightRadii = zeros(size(levels));
relativeDensity = rho / peak;
for j = 1:numel(levels)
    leftIndex = find(relativeDensity(peakIndex:-1:1) < levels(j), 1, 'first');
    if isempty(leftIndex)
        leftRadii(j) = grid.L;
    else
        actualIndex = peakIndex - leftIndex + 1;
        leftRadii(j) = abs(grid.x(actualIndex) - grid.x(peakIndex));
    end
    rightIndex = find(relativeDensity(peakIndex:end) < levels(j), 1, 'first');
    if isempty(rightIndex)
        rightRadii(j) = grid.L;
    else
        actualIndex = peakIndex + rightIndex - 1;
        rightRadii(j) = abs(grid.x(actualIndex) - grid.x(peakIndex));
    end
    radii(j) = max(leftRadii(j), rightRadii(j));
end
diagnostic.relative_levels = levels;
diagnostic.radius_by_level = radii;
diagnostic.radius_left_by_level = leftRadii;
diagnostic.radius_right_by_level = rightRadii;
diagnostic.radius_rel_1e3 = radii(1);
diagnostic.radius_rel_1e6 = radii(2);
diagnostic.radius_rel_1e10 = radii(3);
diagnostic.radius_rel_1e14 = radii(4);
end
