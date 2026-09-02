function [rhoRefN, hatRefN, info] = FourierProjectReference( ...
    rhoRef, gridRef, gridN)
%FOURIERPROJECTREFERENCE Orthogonally restrict a real reference spectrum.
%
% Interior +/- modes are copied exactly. At the even-grid Nyquist mode,
% the +/- reference coefficients are combined into the single real cosine
% degree of freedom represented by the coarse nodal grid. The diagnostic
% projection is never modified to restore positivity.

validateGrids(rhoRef, gridRef, gridN);
[hatRef, modesRef] = src.discretization.ps.FourierCoefficients(rhoRef);
N = gridN.N;
modesN = (-N/2:N/2-1)';
hatRefN = complex(zeros(N, 1));
interior = abs(modesN) < N / 2;
[present, locations] = ismember(modesN(interior), modesRef);
if ~all(present)
    error('src:diagnostics:FourierProjectReference:NonNestedModes', ...
        'The reference grid does not contain every coarse Fourier mode.');
end
hatRefN(interior) = hatRef(locations);
negativeNyquist = hatRef(modesRef == -N / 2);
positiveNyquist = hatRef(modesRef == N / 2);
hatRefN(modesN == -N / 2) = negativeNyquist + positiveNyquist;
rhoRefN = src.discretization.ps.ValuesFromCoefficients(hatRefN);
if ~isreal(rhoRefN)
    error('src:diagnostics:FourierProjectReference:ComplexProjection', ...
        'A real reference state produced a materially complex projection.');
end

negative = rhoRefN < 0;
info.modes = modesN;
info.reference_coefficients = hatRef;
info.reference_modes = modesRef;
info.min_projected_reference = min(rhoRefN);
info.negative_count = nnz(negative);
if any(negative)
    info.negative_min = min(rhoRefN(negative));
else
    info.negative_min = 0;
end
info.nyquist_combined_coefficient = hatRefN(1);
info.grid_ratio = gridRef.N / gridN.N;
end

function validateGrids(rhoRef, gridRef, gridN)
required = {'N', 'L', 'h', 'x'};
if ~all(isfield(gridRef, required)) || ~all(isfield(gridN, required))
    error('src:diagnostics:FourierProjectReference:InvalidGrid', ...
        'Both grids must contain N, L, h, and x.');
end
if numel(rhoRef) ~= gridRef.N || gridRef.N <= gridN.N ...
        || mod(gridRef.N, 2) ~= 0 || mod(gridN.N, 2) ~= 0
    error('src:diagnostics:FourierProjectReference:InvalidSizes', ...
        'Require even grids with numel(rhoRef)=Nref and Nref>N.');
end
scale = max([1, abs(gridRef.L), abs(gridN.L)]);
if abs(gridRef.L - gridN.L) > 100 * eps(scale)
    error('src:diagnostics:FourierProjectReference:DomainMismatch', ...
        'Reference and coarse grids must use the same L.');
end
if mod(gridRef.N, gridN.N) ~= 0
    error('src:diagnostics:FourierProjectReference:NonNestedGrids', ...
        'Strict comparison requires integer Nref/N.');
end
expectedRef = -gridRef.L + (0:gridRef.N-1)' * gridRef.h;
expectedN = -gridN.L + (0:gridN.N-1)' * gridN.h;
if max(abs(gridRef.x(:) - expectedRef)) > 100 * eps(scale) ...
        || max(abs(gridN.x(:) - expectedN)) > 100 * eps(scale)
    error('src:diagnostics:FourierProjectReference:GridConventionMismatch', ...
        'Both grids must use the periodic convention x=-L+(0:N-1)h.');
end
end
