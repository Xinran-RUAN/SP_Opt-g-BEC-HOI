function stats = test_convexity_all_potentials()
%TEST_CONVEXITY_ALL_POTENTIALS Random feasible segment checks.

parameters = model.DefaultParameters1D();
parameters.L = 6;
parameters.N = 64;
grid = model.SetupGrid1D(parameters);
problem.grid = grid;
problem.plan = src.discretization.ps.Plan1D(grid);
problem.V = model.BuildPotential(grid, parameters.potential_label);
problem.beta = parameters.beta;
problem.delta = parameters.delta;
problem.mass = parameters.mass;
problem.regularization.name = 'shift_smooth';
problem.regularization.epsilon = parameters.epsilon;
problem.regularization.transition_width = parameters.epsilon;

rng(97);
names = src.potential.SupportedNames();
maxDefects = zeros(numel(names), 1);
for j = 1:numel(names)
    problem.potential_regularization.name = names{j};
    for trial = 1:40
        rho0 = exp(0.8 * randn(grid.N, 1));
        rho1 = exp(0.8 * randn(grid.N, 1));
        rho0 = rho0 / src.constraints.Mass(rho0, grid.h);
        rho1 = rho1 / src.constraints.Mass(rho1, grid.h);
        theta = rand();
        mixed = theta * rho0 + (1 - theta) * rho1;
        e0 = src.discretization.ps.Energy(rho0, problem);
        e1 = src.discretization.ps.Energy(rho1, problem);
        eMixed = src.discretization.ps.Energy(mixed, problem);
        defect = eMixed - theta * e0 - (1 - theta) * e1;
        maxDefects(j) = max(maxDefects(j), defect);
        tolerance = 2e-11 * max([1, abs(e0), abs(e1)]);
        assert(defect <= tolerance, ...
            'Potential convexity defect for %s is %.3e.', ...
            names{j}, defect);
    end
    fprintf('  potential convexity %-20s max defect %.3e\n', ...
        names{j}, maxDefects(j));
end
stats.names = names;
stats.max_defect = maxDefects;
stats.all_passed = true;
end
