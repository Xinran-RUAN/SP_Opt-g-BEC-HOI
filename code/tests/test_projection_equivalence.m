function stats = test_projection_equivalence()
%TEST_PROJECTION_EQUIVALENCE Compare independent projection backends.

rng(23);
Nlist = [8, 32, 128, 512];
maxDifference = 0;
maxMassError = 0;
minimumDensity = inf;

for N = Nlist
    L = 3.7;
    h = 2 * L / N;
    for sample = 1:12
        z = 0.7 * randn(N, 1) + 0.2 * sin((1:N)' * sample);
        rhoSimplex = src.constraints.ProjectSimplex(z, 1, h);
        rhoSemismooth = src.constraints.ProjectPositiveConservative( ...
            z, 1, h, 1e-14);
        maxDifference = max(maxDifference, norm(rhoSimplex - rhoSemismooth));
        maxMassError = max(maxMassError, max( ...
            abs(h * sum(rhoSimplex) - 1), ...
            abs(h * sum(rhoSemismooth) - 1)));
        minimumDensity = min(minimumDensity, min([rhoSimplex; rhoSemismooth]));
    end
end

assert(maxDifference <= 1e-12, ...
    'Projection difference %.3e exceeds tolerance.', maxDifference);
assert(maxMassError <= 1e-12, ...
    'Projection mass error %.3e exceeds tolerance.', maxMassError);
assert(minimumDensity >= -1e-14, ...
    'Projection produced density %.3e.', minimumDensity);

stats.max_difference = maxDifference;
stats.max_mass_error = maxMassError;
stats.min_density = minimumDensity;
fprintf('test_projection_equivalence: difference %.3e, mass %.3e\n', ...
    maxDifference, maxMassError);
end
