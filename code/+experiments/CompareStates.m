function comparison = CompareStates(rhoCoarse, rhoReference, referenceGrid)
%COMPARESTATES Secondary prolonged-to-reference-grid state diagnostic.
%
% Primary mesh diagnostics must use src.diagnostics.SpectralStateComparison
% (rho_N versus P_N rho_ref). This helper is retained for visualization and
% backward-compatible secondary output only.

rhoOnReference = src.discretization.ps.Prolong( ...
    rhoCoarse, referenceGrid.N);
difference = rhoOnReference(:) - rhoReference(:);
comparison.density_L2_error = sqrt(referenceGrid.h * sum(difference .^ 2));
comparison.density_Linf_error = max(abs(difference));
comparison.prolonged_L2_error = comparison.density_L2_error;
comparison.prolonged_Linf_error = comparison.density_Linf_error;
comparison.rho_on_reference = rhoOnReference;
end
