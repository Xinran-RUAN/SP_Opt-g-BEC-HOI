function stats = test_fourier_prolongation()
%TEST_FOURIER_PROLONGATION Exact interpolation of low-frequency data.

L = 5;
Ncoarse = 32;
xCoarse = -L + (0:Ncoarse-1)' * (2 * L / Ncoarse);
fCoarse = testFunction(xCoarse, L);
targets = [64, 128];
errors = zeros(size(targets));

for j = 1:numel(targets)
    Nfine = targets(j);
    xFine = -L + (0:Nfine-1)' * (2 * L / Nfine);
    fFine = src.discretization.ps.Prolong(fCoarse, Nfine);
    errors(j) = max(abs(fFine - testFunction(xFine, L)));
end

assert(max(errors) <= 2e-12, ...
    'Fourier prolongation error %.3e is too large.', max(errors));
stats.Nfine = targets;
stats.errors = errors;
fprintf('test_fourier_prolongation: max error %.3e\n', max(errors));
end

function values = testFunction(x, L)
values = 0.2 + 0.3 * cos(3 * pi * x / L) ...
    - 0.17 * sin(5 * pi * x / L) + 0.08 * cos(7 * pi * x / L);
end
