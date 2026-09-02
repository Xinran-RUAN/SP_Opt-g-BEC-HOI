function diagnostic = ActualEndpointCompatibility( ...
    rho, problem, endpointData)
%ACTUALENDPOINTCOMPATIBILITY Endpoint jumps of rho and V*dp_sigma(rho).
%
% rho and its derivative are evaluated as the production trigonometric
% interpolant. Hence their values at x=-L and x=L coincide to the periodic
% FFT convention. endpointData supplies only the analytic V and V' values.

required = {'V', 'dV'};
if nargin < 3 || ~isstruct(endpointData) ...
        || ~all(isfield(endpointData, required)) ...
        || numel(endpointData.V) ~= 2 || numel(endpointData.dV) ~= 2
    error('src:diagnostics:ActualEndpointCompatibility:InvalidEndpointData', ...
        'endpointData.V and endpointData.dV must be [left,right] pairs.');
end
rho = rho(:);
drho = src.discretization.ps.FirstDerivative(rho, problem.plan);
rhoEndpoints = [rho(1); rho(1)];
drhoEndpoints = [drho(1); drho(1)];
fisher = src.regularization.ResolveFisher(problem);
[~, dp, d2p, ~] = src.potential.Evaluate( ...
    rhoEndpoints, problem.potential_regularization, fisher.epsilon);
V = endpointData.V(:);
dV = endpointData.dV(:);
gV = V .* dp;
dgV = dV .* dp + V .* d2p .* drhoEndpoints;

diagnostic.rho_endpoint_values = rhoEndpoints;
diagnostic.rho_endpoint_derivatives = drhoEndpoints;
diagnostic.rho_value_jump = rhoEndpoints(2) - rhoEndpoints(1);
diagnostic.rho_derivative_jump = drhoEndpoints(2) - drhoEndpoints(1);
diagnostic.GV_endpoint_values = gV;
diagnostic.GV_endpoint_derivatives = dgV;
diagnostic.GV_value_jump = gV(2) - gV(1);
diagnostic.GV_derivative_jump = dgV(2) - dgV(1);
diagnostic.dp_endpoint = dp;
diagnostic.d2p_endpoint = d2p;
diagnostic.nodal_edge_value_mismatch = rho(end) - rho(1);
diagnostic.nodal_edge_derivative_mismatch = drho(end) - drho(1);
end
