function diagnostic = FourierTailDiagnostics(rho)
%FOURIERTAILDIAGNOSTICS Normalized Fourier coefficients and tail energies.

[coefficients, modes] = src.discretization.ps.FourierCoefficients(rho);
power = abs(coefficients) .^ 2;
totalPower = sum(power);
if totalPower == 0
    quarter = 0;
    third = 0;
else
    N = numel(rho);
    quarter = sum(power(abs(modes) >= N / 4)) / totalPower;
    third = sum(power(abs(modes) >= N / 3)) / totalPower;
end
diagnostic.modes = modes;
diagnostic.rho_hat = coefficients;
diagnostic.rho_hat_abs = abs(coefficients);
diagnostic.tail_ratio_quarter = quarter;
diagnostic.tail_ratio_third = third;
fitMask = abs(modes) >= numel(rho) / 8 ...
    & abs(modes) <= numel(rho) / 3 & diagnostic.rho_hat_abs > 0;
if nnz(fitMask) >= 4
    fit = polyfit(abs(modes(fitMask)), ...
        log(diagnostic.rho_hat_abs(fitMask)), 1);
    diagnostic.fourier_decay_slope = fit(1);
    diagnostic.fourier_decay_fit_count = nnz(fitMask);
else
    diagnostic.fourier_decay_slope = NaN;
    diagnostic.fourier_decay_fit_count = nnz(fitMask);
end
end
