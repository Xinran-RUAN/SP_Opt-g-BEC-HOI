function stats = test_entropy_prox_mass()
%TEST_ENTROPY_PROX_MASS Mass accuracy and explicit underflow classification.

rng(31);
N = 73;
h = 0.2;
mass = 1.3;
z = 0.2 * randn(N, 1);
etaValues = [1e-1, 1e-3, 1e-6, 1e-10];
tauValues = [0.1, 1, 3];
maximumMassError = 0;
underflowCases = 0;
for eta = etaValues
    for tau = tauValues
        [rho, info] = src.entropy.SimplexProx( ...
            z, tau * eta, mass, h, struct());
        maximumMassError = max(maximumMassError, info.mass_error);
        assert(info.mass_error <= 1e-12, ...
            'Entropy prox mass error is %.3e.', info.mass_error);
        if info.underflow_count == 0
            assert(all(rho > 0), 'Non-underflow prox result is not positive.');
        else
            underflowCases = underflowCases + 1;
            assert(nnz(rho == 0) == info.underflow_count, ...
                'Underflow count does not match zero entries.');
        end
    end
end
stats.maximum_mass_error = maximumMassError;
stats.underflow_cases = underflowCases;
fprintf('test_entropy_prox_mass: max mass %.3e, underflow cases %d\n', ...
    maximumMassError, underflowCases);
end
