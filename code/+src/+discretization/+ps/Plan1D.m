function plan = Plan1D(grid)
%PLAN1D Create the FFT differentiation plan for a periodic even grid.
%
% MATLAB uses fft(u)_k = sum_j u_j exp(-2*pi*i*j*k/N). The derivative
% multiplier is i*kappa. For even N the Nyquist multiplier is explicitly
% zero, which makes real differentiation skew-adjoint under h*sum(u.*v).

if ~isstruct(grid) || ~all(isfield(grid, {'N', 'L', 'h'}))
    error('src:discretization:ps:Plan1D:InvalidGrid', ...
        'grid must contain N, L, and h.');
end
N = grid.N;
if N ~= round(N) || mod(N, 2) ~= 0
    error('src:discretization:ps:Plan1D:InvalidN', ...
        'N must be even.');
end

integerModes = [0:(N/2-1), 0, (-N/2+1):-1]';
wavenumbers = (pi / grid.L) * integerModes;
plan.N = N;
plan.L = grid.L;
plan.h = grid.h;
plan.domain_length = 2 * grid.L;
plan.integer_modes = integerModes;
plan.wavenumbers = wavenumbers;
plan.first_derivative_multiplier = 1i * wavenumbers;
plan.coefficient_convention = 'c_k = fft(u)_k / N';
plan.nyquist_derivative_is_zero = true;
end
