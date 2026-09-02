function grid = SetupGrid1D(parameters)
%SETUPGRID1D Construct the periodic collocation grid [-L,L).

if ~isstruct(parameters) || ~all(isfield(parameters, {'L', 'N'}))
    error('model:SetupGrid1D:InvalidInput', ...
        'Input must be a structure containing L and N.');
end
L = parameters.L;
N = parameters.N;
if ~isscalar(L) || ~isfinite(L) || L <= 0
    error('model:SetupGrid1D:InvalidL', 'L must be a positive scalar.');
end
if ~isscalar(N) || N ~= round(N) || N < 4 || mod(N, 2) ~= 0
    error('model:SetupGrid1D:InvalidN', 'N must be an even integer >= 4.');
end

h = 2 * L / N;
grid.L = L;
grid.N = N;
grid.h = h;
grid.domain_length = 2 * L;
grid.x = -L + (0:N-1)' * h;
end
