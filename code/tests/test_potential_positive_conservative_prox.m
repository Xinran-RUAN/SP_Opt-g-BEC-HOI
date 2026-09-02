function stats = test_potential_positive_conservative_prox()
%TEST_POTENTIAL_POSITIVE_CONSERVATIVE_PROX Mass, KKT, and limit tests.

rng(211);
names = src.potential.SupportedNames();
Nlist = [9, 32, 97];
tauList = [1e-3, 0.07, 1];
epsilonList = [1e-2, 1e-3];
options.mass_tol = 1e-13;
options.lambda_max_iter = 50;
options.inner_tol = 1e-13;
options.inner_max_iter = 30;
options.projection_tol = 1e-14;

maxMassError = 0;
minimumDensity = inf;
maxKKTResidual = 0;
maxLinearDifference = 0;
maxUniquenessDifference = 0;
for e = 1:numel(epsilonList)
    epsilon = epsilonList(e);
    for n = 1:numel(Nlist)
        N = Nlist(n);
        L = 3.25;
        h = 2 * L / N;
        mass = 0.7;
        x = -L + (0:N-1)' * h;
        V = 0.1 + 0.5 * x .^ 2;
        z = 0.3 * randn(N, 1) + 0.1 * cos((1:N)' * 0.37);
        for t = 1:numel(tauList)
            tau = tauList(t);
            for p = 1:numel(names)
                [potreg, ~] = src.potential.Validate( ...
                    struct('name', names{p}), epsilon);
                [rho, info] = src.potential.PositiveConservativeProx( ...
                    z, tau, V, potreg, mass, h, options);
                maxMassError = max(maxMassError, info.mass_error);
                minimumDensity = min(minimumDensity, min(rho));
                maxKKTResidual = max(maxKKTResidual, info.kkt_residual);
                assert(info.converged, 'Potential prox did not converge.');
                assert(info.mass_error <= 1e-12, ...
                    'Potential prox mass error %.3e.', info.mass_error);
                assert(min(rho) >= -1e-14, ...
                    'Potential prox density %.3e.', min(rho));
                assert(info.kkt_residual <= 2e-12, ...
                    'Potential prox KKT residual %.3e.', info.kkt_residual);

                if strcmp(names{p}, 'linear')
                    rhoReference = ...
                        src.constraints.ProjectPositiveConservative( ...
                        z - tau * V, mass, h, options.projection_tol);
                    maxLinearDifference = max(maxLinearDifference, ...
                        max(norm(rho - rhoReference), ...
                        max(abs(rho - rhoReference))));
                else
                    lowStart = options;
                    lowStart.lambda_initial = min(z) - 100;
                    highStart = options;
                    highStart.lambda_initial = max(z) - eps(max(z));
                    rhoLow = src.potential.PositiveConservativeProx( ...
                        z, tau, V, potreg, mass, h, lowStart);
                    rhoHigh = src.potential.PositiveConservativeProx( ...
                        z, tau, V, potreg, mass, h, highStart);
                    maxUniquenessDifference = max(maxUniquenessDifference, ...
                        norm(rhoLow - rhoHigh));
                end
            end
        end
    end
end

assert(maxLinearDifference <= 1e-12, ...
    'Linear potential prox difference %.3e.', maxLinearDifference);
assert(maxUniquenessDifference <= 2e-11, ...
    'Potential prox initial-lambda difference %.3e.', ...
    maxUniquenessDifference);

% Regression for the safeguarded node solve: the equation target
% z-lambda must remain fixed when the Newton bracket is updated.  This
% mixes a saturated large-r node with a high-curvature r~sigma node.
sigmaRegression = 2e-6;
regressionPotential.sigma = sigmaRegression;
regressionPotential.p_sigma = @(r) ...
    r .^ 2 ./ (hypot(r, sigmaRegression) + sigmaRegression);
regressionPotential.dp_sigma = @(r) ...
    r ./ hypot(r, sigmaRegression);
regressionPotential.d2p_sigma = @(r) ...
    (sigmaRegression ./ hypot(r, sigmaRegression)) .^ 2 ...
    ./ hypot(r, sigmaRegression);
regressionPotential.prox_type = 'generic_convex';
knownRho = [189; 1.108e-6; 0.2; 0.05];
regressionV = [128; 81; 4; 0.5];
knownLambda = -7.3;
regressionZ = knownRho + regressionV ...
    .* regressionPotential.dp_sigma(knownRho) + knownLambda;
regressionOptions = options;
regressionOptions.mass_tol = 1e-12;
regressionOptions.inner_tol = 1e-14;
regressionOptions.inner_max_iter = 60;
[regressionRho, regressionInfo] = ...
    src.potential.PositiveConservativeProx( ...
    regressionZ, 1, regressionV, regressionPotential, ...
    sum(knownRho), 1, regressionOptions);
regressionStateError = norm(regressionRho - knownRho, inf);
assert(regressionStateError <= 1e-12 ...
    && regressionInfo.mass_error <= 1e-12 ...
    && regressionInfo.kkt_residual <= 2e-12, ...
    ['Safeguarded prox regression failed: state %.3e, mass %.3e, ' ...
    'KKT %.3e.'], regressionStateError, regressionInfo.mass_error, ...
    regressionInfo.kkt_residual);

% The sqrt prox must approach the linear-potential prox as sigma decreases.
N = 48;
L = 4;
h = 2 * L / N;
mass = 1;
x = -L + (0:N-1)' * h;
V = 0.2 + x .^ 2;
z = 0.25 * sin(1.7 * x) + 0.08 * randn(N, 1);
tau = 0.13;
linear.name = 'linear';
rhoLinear = src.potential.PositiveConservativeProx( ...
    z, tau, V, linear, mass, h, options);
sigmaList = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8];
limitErrors = zeros(size(sigmaList));
for j = 1:numel(sigmaList)
    sqrtMap.name = 'sqrt_squared_scale';
    sqrtMap.sigma = sigmaList(j);
    rhoSqrt = src.potential.PositiveConservativeProx( ...
        z, tau, V, sqrtMap, mass, h, options);
    limitErrors(j) = norm(rhoSqrt - rhoLinear);
end
assert(limitErrors(end) < limitErrors(1) ...
        && all(diff(limitErrors) <= 1e-11), ...
    'sigma->0 prox errors do not decrease: %s', mat2str(limitErrors, 3));

stats.max_mass_error = maxMassError;
stats.min_density = minimumDensity;
stats.max_kkt_residual = maxKKTResidual;
stats.max_linear_difference = maxLinearDifference;
stats.max_uniqueness_difference = maxUniquenessDifference;
stats.safeguarded_regression_state_error = regressionStateError;
stats.sigma_limit_errors = limitErrors;
fprintf(['test_potential_positive_conservative_prox: mass %.3e, ' ...
    'min %.3e, KKT %.3e, linear diff %.3e, uniqueness %.3e\n'], ...
    maxMassError, minimumDensity, maxKKTResidual, ...
    maxLinearDifference, maxUniquenessDifference);
fprintf('  sigma->0 errors: %s\n', sprintf('%.3e ', limitErrors));
end
