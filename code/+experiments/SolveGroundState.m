function result = SolveGroundState(config, rho0)
%SOLVEGROUNDSTATE Assemble and solve one configured discrete problem.

parameters = config.parameters;
grid = model.SetupGrid1D(parameters);
plan = src.discretization.ps.Plan1D(grid);
if isfield(config, 'trapping_potential') ...
        && ~isempty(config.trapping_potential)
    [V, trappingPotential] = model.EvaluatePotential( ...
        grid, config.trapping_potential);
else
    % Backward-compatible path for archived/legacy configurations.
    V = model.BuildPotential(grid, parameters.potential_label);
    legacyLabel = parameters.potential_label;
    trappingPotential.V = @(x) legacyPotential(x, legacyLabel);
    trappingPotential.label = legacyLabel;
    trappingPotential.values = V;
    trappingPotential.minimum = min(V);
    trappingPotential.maximum = max(V);
    trappingPotential.source = 'legacy_label';
end

if isfield(config, 'fisher_regularization') ...
        && ~isempty(config.fisher_regularization)
    fisherRegularization = config.fisher_regularization;
    regularization.name = 'inline_handles';
    regularization.epsilon = fisherRegularization.epsilon;
    regularization.transition_width = fisherRegularization.epsilon;
    legacyRegularizationInfo.name = 'inline_handles';
    legacyRegularizationInfo.epsilon = fisherRegularization.epsilon;
    legacyRegularizationInfo.is_convex_preserving = true;
else
    regularization = config.regularization;
    regularization.epsilon = parameters.epsilon;
    [regularization, legacyRegularizationInfo] = ...
        src.regularization.Validate(regularization);
    fisherRegularization = src.regularization.MakeBuiltIn(regularization);
end
rhoTest = [0; parameters.mass / grid.domain_length; ...
    parameters.mass / grid.h];
[fisherRegularization, fisherRegularizationInfo] = ...
    src.regularization.ValidateFisher(fisherRegularization, rhoTest);

problem.grid = grid;
problem.plan = plan;
problem.V = V;
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.regularization = regularization;
problem.fisher_regularization = fisherRegularization;
problem.trapping_potential = trappingPotential;
problem.potential_label = trappingPotential.label; % compatibility metadata
if isfield(config, 'potential_regularization')
    problem.potential_regularization = config.potential_regularization;
else
    problem.potential_regularization.name = 'linear';
end
if isfield(config, 'entropy')
    problem.entropy = config.entropy;
else
    problem.entropy.enabled = false;
    problem.entropy.eta = 0;
end
entropyActive = problem.entropy.enabled && problem.entropy.eta > 0;
if entropyActive && strcmpi(config.solver.name, 'spg')
    config.solver.name = 'fista_cd';
end

if nargin < 2 || isempty(rho0)
    rho0 = model.InitialDensity( ...
        grid, parameters.mass, config.solver.projection_name);
end

result = src.SolveGroundState1D(problem, rho0, config.solver);
result.parameters = parameters;
result.regularization = legacyRegularizationInfo;
result.fisher_regularization = fisherRegularizationInfo;
result.grid = grid;
result.potential = V;
result.trapping_potential = trappingPotential;
result.target_energy = result.energy;
result.baseline_energy = src.discretization.ps.BaselineEnergy( ...
    result.rho, result.problem);

tailWidth = min(2, 0.1 * grid.L);
tailMask = abs(grid.x) >= grid.L - tailWidth;
tailMass = grid.h * sum(result.rho(tailMask));
tailMax = max(result.rho(tailMask));
result.diagnostics.tail_width = tailWidth;
result.diagnostics.tail_mass = tailMass;
result.diagnostics.tail_max = tailMax;
result.diagnostics.tail_indicator = max(tailMass, tailMax);
result.diagnostics.rho_left_boundary = result.rho(1);
result.diagnostics.rho_right_boundary = result.rho(end);
result.diagnostics.boundary_value_mismatch = ...
    abs(result.rho(1) - result.rho(end));
derivative = src.discretization.ps.FirstDerivative(result.rho, problem.plan);
result.diagnostics.boundary_derivative_mismatch = ...
    abs(derivative(1) - derivative(end));
result.diagnostics.active_set = ...
    src.diagnostics.ActiveSetDiagnostics(result.rho, grid);
result.diagnostics.fourier_tail = ...
    src.diagnostics.FourierTailDiagnostics(result.rho);
result.diagnostics.far_field = src.diagnostics.FarFieldDiagnostics( ...
    result.rho, grid, problem.plan);
result.diagnostics.relative_support = ...
    src.diagnostics.RelativeSupportDiagnostics(result.rho, grid);
result.diagnostics.trapping_potential_modification = ...
    potentialModificationDiagnostics( ...
    result.rho, grid, V, trappingPotential);
end

function values = legacyPotential(x, label)
legacyGrid.x = x(:);
values = model.BuildPotential(legacyGrid, label);
end

function diagnostic = potentialModificationDiagnostics( ...
    rho, grid, V, potential)
diagnostic.enabled = isfield(potential, 'boundary_periodicized') ...
    && logical(potential.boundary_periodicized);
diagnostic.rho_max = NaN;
diagnostic.mass = NaN;
diagnostic.potential_change_max = NaN;
diagnostic.core_potential_error = NaN;
if ~diagnostic.enabled
    return;
end
required = {'modification_start', 'reference_V'};
if ~all(isfield(potential, required)) ...
        || ~isa(potential.reference_V, 'function_handle')
    error('experiments:SolveGroundState:IncompletePeriodicPotential', ...
        ['A boundary-periodicized potential requires modification_start ' ...
        'and reference_V metadata.']);
end
R0 = potential.modification_start;
modified = abs(grid.x) >= R0;
core = abs(grid.x) <= R0;
referenceValues = potential.reference_V(grid.x);
referenceValues = referenceValues(:);
diagnostic.modification_start = R0;
diagnostic.transition_end = potential.transition_end;
diagnostic.rho_max = max(rho(modified));
diagnostic.mass = grid.h * sum(rho(modified));
diagnostic.potential_change_max = max(abs(V(modified) ...
    - referenceValues(modified)));
diagnostic.core_potential_error = max(abs(V(core) ...
    - referenceValues(core)));
if isfield(potential, 'negligible_density_tol')
    diagnostic.density_tol = potential.negligible_density_tol;
else
    diagnostic.density_tol = NaN;
end
if isfield(potential, 'negligible_mass_tol')
    diagnostic.mass_tol = potential.negligible_mass_tol;
else
    diagnostic.mass_tol = NaN;
end
diagnostic.density_negligible = isnan(diagnostic.density_tol) ...
    || diagnostic.rho_max <= diagnostic.density_tol;
diagnostic.mass_negligible = isnan(diagnostic.mass_tol) ...
    || diagnostic.mass <= diagnostic.mass_tol;
end
