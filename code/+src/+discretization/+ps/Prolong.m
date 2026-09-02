function rhoFine = Prolong(rhoCoarse, Nfine)
%PROLONG Fourier zero-padding interpolation to a finer even grid.
%
% Coefficients use c_k = fft(rho)_k/N. For an even coarse grid the real
% Nyquist coefficient is split equally between the +/- Nyquist modes on
% the fine grid. This gives the real trigonometric interpolant.

wasRow = isrow(rhoCoarse);
rhoCoarse = rhoCoarse(:);
N = numel(rhoCoarse);
if mod(N, 2) ~= 0 || ~isscalar(Nfine) || Nfine ~= round(Nfine) || ...
        mod(Nfine, 2) ~= 0 || Nfine < N
    error('src:discretization:ps:Prolong:InvalidSize', ...
        'Coarse and fine sizes must be even integers with Nfine >= N.');
end
if Nfine == N
    rhoFine = rhoCoarse;
else
    [coefficients, modes] = ...
        src.discretization.ps.FourierCoefficients(rhoCoarse);
    modesFine = (-Nfine/2:Nfine/2-1)';
    coefficientsFine = complex(zeros(Nfine, 1));
    interior = abs(modes) < N / 2;
    [present, locations] = ismember(modes(interior), modesFine);
    if ~all(present)
        error('src:discretization:ps:Prolong:ModeMappingFailure', ...
            'Coarse interior modes are absent from the fine grid.');
    end
    coefficientsFine(locations) = coefficients(interior);
    nyquist = coefficients(modes == -N / 2);
    coefficientsFine(modesFine == -N / 2) = nyquist / 2;
    coefficientsFine(modesFine == N / 2) = nyquist / 2;
    rhoFine = src.discretization.ps.ValuesFromCoefficients(coefficientsFine);
end
if wasRow
    rhoFine = rhoFine.';
end
end
