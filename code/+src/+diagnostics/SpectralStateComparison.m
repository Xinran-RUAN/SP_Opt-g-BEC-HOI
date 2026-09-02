function comparison = SpectralStateComparison(rhoN, gridN, rhoRef, gridRef)
%SPECTRALSTATECOMPARISON Primary rho_N versus P_N rho_ref comparison.

rhoN = rhoN(:);
if numel(rhoN) ~= gridN.N
    error('src:diagnostics:SpectralStateComparison:SizeMismatch', ...
        'rhoN must contain gridN.N values.');
end
[rhoRefN, hatRefN, projectionInfo] = ...
    src.diagnostics.FourierProjectReference(rhoRef, gridRef, gridN);
difference = rhoN - rhoRefN;
comparison.resolved_L2_error = sqrt(gridN.h * sum(abs(difference) .^ 2));
comparison.resolved_Linf_error = max(abs(difference));

[hatN, modesN] = src.discretization.ps.FourierCoefficients(rhoN);
if ~isequal(modesN, projectionInfo.modes)
    error('src:diagnostics:SpectralStateComparison:ModeMismatch', ...
        'The coefficient and projection mode conventions disagree.');
end
coefficientDifference = hatN - hatRefN;
comparison.fourier_resolved_coeff_error = ...
    sqrt(sum(abs(coefficientDifference) .^ 2));
comparison.parseval_resolved_L2 = sqrt(gridN.domain_length) ...
    * comparison.fourier_resolved_coeff_error;
comparison.parseval_consistency_error = abs( ...
    comparison.resolved_L2_error - comparison.parseval_resolved_L2);

tail = src.diagnostics.FourierTailError(rhoRef, gridRef, gridN);
comparison.reference_tail_L2 = tail.reference_tail_L2;
comparison.reference_tail_energy_fraction = ...
    tail.reference_tail_energy_fraction;
comparison.total_spectral_L2_error = hypot( ...
    comparison.resolved_L2_error, comparison.reference_tail_L2);

rhoNOnRef = src.discretization.ps.Prolong(rhoN, gridRef.N);
prolongedDifference = rhoNOnRef - rhoRef(:);
comparison.prolonged_L2_error = sqrt( ...
    gridRef.h * sum(abs(prolongedDifference) .^ 2));
comparison.prolonged_Linf_error = max(abs(prolongedDifference));
comparison.rho_reference_projected_to_N = rhoRefN;
comparison.reference_projection_info = projectionInfo;
comparison.rhoN_on_reference_grid = rhoNOnRef;

if comparison.reference_tail_L2 == 0
    ratio = Inf;
else
    ratio = comparison.resolved_L2_error / comparison.reference_tail_L2;
end
comparison.resolved_to_tail_ratio = ratio;
if ratio > 10
    comparison.regime = 'resolved_error_dominated';
elseif ratio < 0.1
    comparison.regime = 'reference_tail_dominated';
else
    comparison.regime = 'balanced';
end
end
