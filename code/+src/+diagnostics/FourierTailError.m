function diagnostic = FourierTailError(rhoRef, gridRef, gridN)
%FOURIERTAILERROR Reference-space omitted-mode error from Parseval.

[~, hatRefN, projection] = ...
    src.diagnostics.FourierProjectReference(rhoRef, gridRef, gridN);
hatRef = projection.reference_coefficients;
modesRef = projection.reference_modes;
N = gridN.N;
modesN = projection.modes;

hatProjectedOnRef = complex(zeros(gridRef.N, 1));
interior = abs(modesN) < N / 2;
[~, locations] = ismember(modesN(interior), modesRef);
hatProjectedOnRef(locations) = hatRefN(interior);
nyquist = hatRefN(modesN == -N / 2);
hatProjectedOnRef(modesRef == -N / 2) = nyquist / 2;
hatProjectedOnRef(modesRef == N / 2) = nyquist / 2;

tailCoefficients = hatRef - hatProjectedOnRef;
tailPower = sum(abs(tailCoefficients) .^ 2);
totalPower = sum(abs(hatRef) .^ 2);
diagnostic.reference_tail_L2 = sqrt(gridRef.domain_length * tailPower);
if totalPower == 0
    diagnostic.reference_tail_energy_fraction = 0;
else
    diagnostic.reference_tail_energy_fraction = tailPower / totalPower;
end
diagnostic.tail_coefficients = tailCoefficients;
diagnostic.projected_coefficients_on_reference_grid = hatProjectedOnRef;
end
