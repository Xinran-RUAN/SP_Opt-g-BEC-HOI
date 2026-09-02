function diagnostic = ActiveSetDiagnostics(rho, grid, options)
%ACTIVESETDIAGNOSTICS Threshold-resolved zero and support diagnostics.

if nargin < 3
    options = struct();
end
if ~isfield(options, 'zero_tol_list') || isempty(options.zero_tol_list)
    options.zero_tol_list = [0, 1e-15, 1e-12, 1e-10, 1e-8];
end
rho = rho(:);
tolerances = options.zero_tol_list(:).';
counts = zeros(size(tolerances));
fractions = zeros(size(tolerances));
for j = 1:numel(tolerances)
    counts(j) = nnz(rho <= tolerances(j));
    fractions(j) = counts(j) / numel(rho);
end
positive = rho > 0;
diagnostic.zero_tol_list = tolerances;
diagnostic.zero_count_by_tol = counts;
diagnostic.zero_fraction_by_tol = fractions;
diagnostic.exact_zero_count = counts(1);
[~, numericalIndex] = min(abs(tolerances - 1e-12));
diagnostic.numerical_zero_count = counts(numericalIndex);
diagnostic.active_fraction = fractions(numericalIndex);
diagnostic.min_density = min(rho);
if any(positive)
    diagnostic.min_positive_density = min(rho(positive));
else
    diagnostic.min_positive_density = NaN;
end

if ~isfield(options, 'support_tol') || isempty(options.support_tol)
    options.support_tol = 1e-12 * max(rho);
end
support = find(rho > options.support_tol);
diagnostic.support_tol = options.support_tol;
if isempty(support)
    diagnostic.support_left = NaN;
    diagnostic.support_right = NaN;
    diagnostic.support_left_x = NaN;
    diagnostic.support_right_x = NaN;
    diagnostic.support_width = 0;
    diagnostic.support_left_fourier_derivative = NaN;
    diagnostic.support_right_fourier_derivative = NaN;
    diagnostic.support_left_fd_slope_jump = NaN;
    diagnostic.support_right_fd_slope_jump = NaN;
    return;
end

left = support(1);
right = support(end);
diagnostic.support_left = left;
diagnostic.support_right = right;
diagnostic.support_left_x = grid.x(left);
diagnostic.support_right_x = grid.x(right);
diagnostic.support_width = grid.x(right) - grid.x(left);
plan = src.discretization.ps.Plan1D(grid);
derivative = src.discretization.ps.FirstDerivative(rho, plan);
diagnostic.support_left_fourier_derivative = derivative(left);
diagnostic.support_right_fourier_derivative = derivative(right);
diagnostic.support_left_fd_slope_jump = localSlopeJump(rho, left, grid.h);
diagnostic.support_right_fd_slope_jump = localSlopeJump(rho, right, grid.h);
end

function value = localSlopeJump(rho, index, h)
N = numel(rho);
previous = mod(index - 2, N) + 1;
next = mod(index, N) + 1;
value = abs((rho(next) - rho(index)) / h ...
    - (rho(index) - rho(previous)) / h);
end
