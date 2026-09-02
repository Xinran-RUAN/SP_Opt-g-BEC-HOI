function stats = test_harmonic_cinf_periodic_potential()
%TEST_HARMONIC_CINF_PERIODIC_POTENTIAL Core identity and flat boundary.

L = 8;
R0 = 0.75 * L;
R1 = 0.90 * L;
x = linspace(-L, L, 4001)';
[V, info] = model.HarmonicCInfPeriodicPotential(x, L, R0, R1);
harmonic = 0.5 * x .^ 2;
core = abs(x) <= R0;
flat = abs(x) >= R1;
positiveSide = x >= 0;

coreError = max(abs(V(core) - harmonic(core)));
flatError = max(abs(V(flat) - 0.5 * L ^ 2));
symmetryError = max(abs(V - flipud(V)));
minimumIncrement = min(diff(V(positiveSide)));
valueJump = abs(V(1) - V(end));
dx = x(2) - x(1);
leftDerivative = (V(2) - V(1)) / dx;
rightDerivative = (V(end) - V(end-1)) / dx;
derivativeJump = abs(rightDerivative - leftDerivative);

assert(coreError <= 10 * eps(max(1, max(harmonic))), ...
    'Harmonic core was modified by %.3e.', coreError);
assert(flatError == 0, 'Boundary plateau error %.3e.', flatError);
assert(info.minimum >= 0 && minimumIncrement >= -1e-12, ...
    'Continuation created a negative value or an artificial boundary well.');
assert(symmetryError <= 100 * eps(max(V)) ...
    && valueJump == 0 && derivativeJump == 0, ...
    'Periodic endpoint compatibility failed.');

stats.core_error = coreError;
stats.flat_error = flatError;
stats.symmetry_error = symmetryError;
stats.value_jump = valueJump;
stats.derivative_jump = derivativeJump;
stats.minimum_increment = minimumIncrement;
fprintf(['test_harmonic_cinf_periodic_potential: core %.3e, flat %.3e, ' ...
    'value/derivative jump %.3e/%.3e\n'], ...
    coreError, flatError, valueJump, derivativeJump);
end
