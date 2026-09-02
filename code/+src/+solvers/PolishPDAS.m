function result = PolishPDAS(rho0, problem, options)
%POLISHPDAS Compatibility entry that explicitly selects PDAS-GMRES.

options.linear_solver = 'pdas_gmres';
result = src.solvers.PolishKKT(rho0, problem, options);
end
