function stats = test_inline_trapping_potential_interface()
%TEST_INLINE_TRAPPING_POTENTIAL_INTERFACE V(x) comes from run-level handles.

parameters = model.DefaultParameters1D();
parameters.L = 3;
parameters.N = 32;
grid = model.SetupGrid1D(parameters);

potential.V = @(x) 1 + 0.2 * cos(pi * x / parameters.L);
potential.label = 'V(x) = 1 + 0.2 cos(pi x/L)';
[values, normalized] = model.EvaluatePotential(grid, potential);
expected = potential.V(grid.x);
evaluationError = max(abs(values - expected));
assert(evaluationError <= 10 * eps, ...
    'Inline trapping-potential evaluation error %.3e.', evaluationError);
assert(strcmp(normalized.label, potential.label) ...
    && isa(normalized.V, 'function_handle'));

% Verify that the configured expression, rather than the legacy harmonic
% label, is assembled into problem.V.
config = experiments.DefaultConfig();
config.parameters.L = parameters.L;
config.parameters.N = parameters.N;
config.trapping_potential = potential;
config.solver.name = 'spg';
config.solver.max_iter = 2;
config.solver.polish_mode = 'none';
config.solver.display = false;
result = experiments.SolveGroundState(config);
assemblyError = max(abs(result.problem.V - expected));
assert(assemblyError <= 10 * eps, ...
    'SolveGroundState did not use inline V(x): %.3e.', assemblyError);
assert(strcmp(result.trapping_potential.label, potential.label));

% Common-grid diagnostics must reevaluate exactly the same V handle.
rho = ones(grid.N, 1) / grid.domain_length;
Ndiag = 2 * grid.N;
common = src.diagnostics.CommonGridEnergy( ...
    rho, grid, result.problem, Ndiag);
diagParameters.L = grid.L;
diagParameters.N = Ndiag;
gridDiag = model.SetupGrid1D(diagParameters);
problemDiag = result.problem;
problemDiag.grid = gridDiag;
problemDiag.plan = src.discretization.ps.Plan1D(gridDiag);
[problemDiag.V, problemDiag.trapping_potential] = ...
    model.EvaluatePotential(gridDiag, potential);
rhoDiag = src.discretization.ps.Prolong(rho, Ndiag);
expectedCommonEnergy = src.discretization.ps.Energy(rhoDiag, problemDiag);
commonGridError = abs(common.common_physical_energy ...
    - expectedCommonEnergy);
assert(common.common_physical_energy_valid && commonGridError <= 1e-13, ...
    'Common-grid energy did not reuse inline V(x): %.3e.', ...
    commonGridError);

stats.evaluation_error = evaluationError;
stats.assembly_error = assemblyError;
stats.common_grid_error = commonGridError;
fprintf(['test_inline_trapping_potential_interface: evaluation %.3e, ' ...
    'assembly %.3e, common-grid %.3e\n'], ...
    evaluationError, assemblyError, commonGridError);
end
