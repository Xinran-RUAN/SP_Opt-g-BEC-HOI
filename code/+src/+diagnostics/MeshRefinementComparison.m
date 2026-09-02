function record = MeshRefinementComparison(resultN, resultRef, Ndiag, commonRef)
%MESHREFINEMENTCOMPARISON Unified state and energy refinement diagnostics.

if nargin < 4 || isempty(commonRef)
    commonRef = src.diagnostics.CommonGridEnergy( ...
        resultRef.rho, resultRef.grid, resultRef.problem, Ndiag);
end
state = src.diagnostics.SpectralStateComparison( ...
    resultN.rho, resultN.grid, resultRef.rho, resultRef.grid);
commonN = src.diagnostics.CommonGridEnergy( ...
    resultN.rho, resultN.grid, resultN.problem, Ndiag);

record.N = resultN.grid.N;
record.native_augmented_energy_error = abs( ...
    resultN.augmented_energy - resultRef.augmented_energy);
record.native_physical_energy_difference = ...
    resultN.physical_energy - resultRef.physical_energy;
record.common_physical_energy = commonN.common_physical_energy;
record.common_physical_energy_valid = ...
    commonN.common_physical_energy_valid ...
    && commonRef.common_physical_energy_valid;
if record.common_physical_energy_valid
    record.common_physical_energy_error = abs( ...
        commonN.common_physical_energy ...
        - commonRef.common_physical_energy);
else
    record.common_physical_energy_error = NaN;
end
record.common_augmented_energy = commonN.common_augmented_energy;
record.common_augmented_energy_valid = ...
    commonN.common_augmented_energy_valid ...
    && commonRef.common_augmented_energy_valid;
if record.common_augmented_energy_valid
    record.common_augmented_energy_error = abs( ...
        commonN.common_augmented_energy ...
        - commonRef.common_augmented_energy);
else
    record.common_augmented_energy_error = NaN;
end
record.min_rho_diag = commonN.min_rho_diag;
record.negative_count_diag = commonN.negative_count_diag;
record.zero_count_diag = commonN.zero_count_diag;

record.resolved_L2_error = state.resolved_L2_error;
record.resolved_Linf_error = state.resolved_Linf_error;
record.reference_tail_L2 = state.reference_tail_L2;
record.reference_tail_energy_fraction = ...
    state.reference_tail_energy_fraction;
record.total_spectral_L2_error = state.total_spectral_L2_error;
record.fourier_resolved_coeff_error = ...
    state.fourier_resolved_coeff_error;
record.parseval_consistency_error = state.parseval_consistency_error;
record.prolonged_L2_error = state.prolonged_L2_error;
record.prolonged_Linf_error = state.prolonged_Linf_error;
record.resolved_to_tail_ratio = state.resolved_to_tail_ratio;
record.regime = state.regime;
record.projected_reference_min = ...
    state.reference_projection_info.min_projected_reference;
record.projected_reference_negative_count = ...
    state.reference_projection_info.negative_count;
record.projected_reference_negative_min = ...
    state.reference_projection_info.negative_min;

record.pg_residual = resultN.diagnostics.pg_residual;
record.kkt_residual = resultN.diagnostics.kkt_residual;
record.fft_tail_ratio = ...
    resultN.diagnostics.fourier_tail.tail_ratio_quarter;
record.tail_mass = resultN.diagnostics.tail_mass;
record.tail_max = resultN.diagnostics.tail_max;
end
