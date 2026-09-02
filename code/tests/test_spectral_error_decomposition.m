function stats = test_spectral_error_decomposition()
%TEST_SPECTRAL_ERROR_DECOMPOSITION Parseval resolved/tail orthogonality.

rng(83);
L = pi;
N = 32;
Nref = 128;
gridN = makeGrid(L, N);
gridRef = makeGrid(L, Nref);
modesN = (-N/2:N/2-1)';
modesRef = (-Nref/2:Nref/2-1)';

hatN = randomRealSpectrum(modesN, 12);
hatN(modesN == -N / 2) = 0;
hatRef = randomRealSpectrum(modesRef, 50);
hatRef(abs(modesRef) == N / 2) = 0;
rhoN = src.discretization.ps.ValuesFromCoefficients(hatN);
rhoRef = src.discretization.ps.ValuesFromCoefficients(hatRef);

comparison = src.diagnostics.SpectralStateComparison( ...
    rhoN, gridN, rhoRef, gridRef);
rhoNOnRef = src.discretization.ps.Prolong(rhoN, Nref);
directSquared = gridRef.h * sum(abs(rhoNOnRef - rhoRef) .^ 2);
decomposedSquared = comparison.resolved_L2_error ^ 2 ...
    + comparison.reference_tail_L2 ^ 2;
scale = max([directSquared, decomposedSquared, eps]);
stats.decomposition_relative_error = ...
    abs(directSquared - decomposedSquared) / scale;
stats.resolved_parseval_relative_error = ...
    comparison.parseval_consistency_error ...
    / max(comparison.resolved_L2_error, eps);
stats.maximum_relative_error = max( ...
    stats.decomposition_relative_error, ...
    stats.resolved_parseval_relative_error);
assert(stats.maximum_relative_error <= 2e-13, ...
    'Spectral Parseval decomposition error is %.3e.', ...
    stats.maximum_relative_error);
fprintf('test_spectral_error_decomposition: max relative %.3e\n', ...
    stats.maximum_relative_error);
end

function coefficients = randomRealSpectrum(modes, maximumMode)
coefficients = complex(zeros(size(modes)));
coefficients(modes == 0) = randn();
for k = 1:maximumMode
    if any(modes == k) && any(modes == -k)
        value = randn() + 1i * randn();
        coefficients(modes == k) = value;
        coefficients(modes == -k) = conj(value);
    end
end
nyquist = min(modes);
coefficients(modes == nyquist) = randn();
end

function grid = makeGrid(L, N)
parameters.L = L;
parameters.N = N;
grid = model.SetupGrid1D(parameters);
end
