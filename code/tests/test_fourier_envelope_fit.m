function stats = test_fourier_envelope_fit()
%TEST_FOURIER_ENVELOPE_FIT Recover known descriptive coefficient laws.

N = 1024;
modes = (-N/2:N/2-1)';
algebraicCoefficients = zeros(N, 1);
nonzero = modes ~= 0;
algebraicCoefficients(nonzero) = 1 ./ abs(modes(nonzero)) .^ 2;
algebraicCoefficients(modes == 0) = 1;
algebraicValues = src.discretization.ps.ValuesFromCoefficients( ...
    algebraicCoefficients);
options.fit_range = [8, 256];
options.relative_floor = 1e-14;
algebraic = src.diagnostics.FourierEnvelopeFit( ...
    algebraicValues, options);

exponentialCoefficients = exp(-0.05 * abs(modes));
exponentialValues = src.discretization.ps.ValuesFromCoefficients( ...
    exponentialCoefficients);
exponential = src.diagnostics.FourierEnvelopeFit( ...
    exponentialValues, options);

stats.algebraic_slope_error = abs(algebraic.algebraic_slope - 2);
stats.algebraic_R2 = algebraic.R2_algebraic;
stats.exponential_slope_error = abs(exponential.exponential_slope - 0.05);
stats.exponential_R2 = exponential.R2_exponential;
assert(stats.algebraic_slope_error <= 1e-10);
assert(stats.algebraic_R2 >= 1 - 1e-12);
assert(stats.exponential_slope_error <= 1e-10);
assert(stats.exponential_R2 >= 1 - 1e-12);
fprintf(['test_fourier_envelope_fit: algebraic m error %.3e, ' ...
    'exponential c error %.3e\n'], ...
    stats.algebraic_slope_error, stats.exponential_slope_error);
end
