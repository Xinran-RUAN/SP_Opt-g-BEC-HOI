function stats = test_spectral_reference_projection()
%TEST_SPECTRAL_REFERENCE_PROJECTION Exact restriction without aliasing.

L = pi;
gridRef = makeGrid(L, 512);
xRef = gridRef.x;
lowReference = 1 + 0.1 * cos(2 * xRef) + 0.03 * sin(5 * xRef);
Nvalues = [32, 64, 128];
lowErrors = zeros(size(Nvalues));
highErrors = zeros(size(Nvalues));
highReference = cos(3 * xRef) + 0.2 * cos(70 * xRef);

for j = 1:numel(Nvalues)
    gridN = makeGrid(L, Nvalues(j));
    lowProjected = src.diagnostics.FourierProjectReference( ...
        lowReference, gridRef, gridN);
    lowExact = 1 + 0.1 * cos(2 * gridN.x) + 0.03 * sin(5 * gridN.x);
    lowErrors(j) = max(abs(lowProjected - lowExact));

    highProjected = src.diagnostics.FourierProjectReference( ...
        highReference, gridRef, gridN);
    highExact = cos(3 * gridN.x);
    highErrors(j) = max(abs(highProjected - highExact));
end
stats.maximum_low_mode_error = max(lowErrors);
stats.maximum_no_alias_error = max(highErrors);
stats.maximum_error = max([lowErrors, highErrors]);
assert(stats.maximum_error <= 5e-14, ...
    'Spectral reference projection error is %.3e.', stats.maximum_error);
fprintf('test_spectral_reference_projection: max error %.3e\n', ...
    stats.maximum_error);
end

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end
