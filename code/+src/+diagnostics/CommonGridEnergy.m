function diagnostic = CommonGridEnergy(rho, sourceGrid, problem, Ndiag)
%COMMONGRIDENERGY Evaluate an unchanged trigonometric state with guards.

if mod(Ndiag, 2) ~= 0 || Ndiag < sourceGrid.N ...
        || mod(Ndiag, sourceGrid.N) ~= 0
    error('src:diagnostics:CommonGridEnergy:InvalidDiagnosticGrid', ...
        'Ndiag must be an even integer multiple of the source N.');
end
rhoDiag = src.discretization.ps.Prolong(rho, Ndiag);
parameters.L = sourceGrid.L;
parameters.N = Ndiag;
gridDiag = model.SetupGrid1D(parameters);
problemDiag = problem;
problemDiag.grid = gridDiag;
problemDiag.plan = src.discretization.ps.Plan1D(gridDiag);
if isfield(problem, 'trapping_potential') ...
        && isfield(problem.trapping_potential, 'V') ...
        && isa(problem.trapping_potential.V, 'function_handle')
    [problemDiag.V, problemDiag.trapping_potential] = ...
        model.EvaluatePotential(gridDiag, problem.trapping_potential);
elseif isfield(problem, 'potential_label')
    problemDiag.V = model.BuildPotential( ...
        gridDiag, problem.potential_label);
else
    problemDiag.V = model.BuildPotential(gridDiag, 'harmonic');
end

diagnostic.N_diag = Ndiag;
diagnostic.rho_diag = rhoDiag;
diagnostic.min_rho_diag = min(rhoDiag);
diagnostic.negative_count_diag = nnz(rhoDiag < 0);
diagnostic.zero_count_diag = nnz(rhoDiag == 0);
diagnostic.common_physical_energy = NaN;
diagnostic.common_physical_energy_valid = false;
diagnostic.common_augmented_energy = NaN;
diagnostic.common_augmented_energy_valid = false;

fisher = src.regularization.ResolveFisher(problem);
rDiag = fisher.r_epsilon(rhoDiag);
physicalValid = all(isfinite(rDiag)) && all(rDiag > 0);
if isfield(fisher, 'legacy_name') ...
        && strcmpi(fisher.legacy_name, 'shift_smooth')
    physicalValid = diagnostic.min_rho_diag ...
        > -fisher.epsilon;
end
if physicalValid
    diagnostic.common_physical_energy = ...
        commonTargetEnergy(rhoDiag, problemDiag);
    diagnostic.common_physical_energy_valid = true;
end

entropyActive = isfield(problem, 'entropy') ...
    && problem.entropy.enabled && problem.entropy.eta > 0;
if entropyActive
    entropyAdmissible = all(rhoDiag > 0);
    if physicalValid && entropyAdmissible
        entropyValue = src.entropy.Value(rhoDiag, gridDiag);
        diagnostic.common_augmented_energy = ...
            diagnostic.common_physical_energy ...
            + problem.entropy.eta * entropyValue;
        diagnostic.common_augmented_energy_valid = true;
    end
elseif physicalValid
    diagnostic.common_augmented_energy = diagnostic.common_physical_energy;
    diagnostic.common_augmented_energy_valid = true;
end
end

function energy = commonTargetEnergy(rho, problem)
q = src.discretization.ps.FirstDerivative(rho, problem.plan);
fisher = src.regularization.ResolveFisher(problem);
r = fisher.r_epsilon(rho);
[potentialDensity, ~] = src.potential.Evaluate( ...
    rho, potentialRegularization(problem), fisher.epsilon);
energy = problem.grid.h * sum(q .^ 2 ./ (8 * r)) ...
    + problem.grid.h * sum(problem.V(:) .* potentialDensity) ...
    + problem.grid.h * (problem.beta / 2) * sum(rho .^ 2) ...
    + problem.grid.h * (problem.delta / 2) * sum(q .^ 2);
end

function potreg = potentialRegularization(problem)
if isfield(problem, 'potential_regularization')
    potreg = problem.potential_regularization;
else
    potreg = src.potential.MakeLinear();
end
end
