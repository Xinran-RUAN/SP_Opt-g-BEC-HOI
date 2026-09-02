function [tf, diagnostics] = IsFeasible(rho, mass, h, tolerance)
%ISFEASIBLE Check nodal positivity and discrete mass conservation.

if nargin < 4 || isempty(tolerance)
    tolerance = 1e-12;
end
diagnostics.mass_error = abs(src.constraints.Mass(rho, h) - mass);
diagnostics.min_density = min(rho(:));
tf = diagnostics.mass_error <= tolerance && ...
    diagnostics.min_density >= -tolerance;
end
