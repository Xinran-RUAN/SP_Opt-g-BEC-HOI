function stats = test_potential_regularization_derivative()
%TEST_POTENTIAL_REGULARIZATION_DERIVATIVE Centered checks of p'(rho).

epsilon = 1e-3;
rho = [0; logspace(-9, 0, 80)'];
names = src.potential.SupportedNames();
bestErrors = zeros(numel(names), 1);
for j = 1:numel(names)
    potreg.name = names{j};
    [~, derivative, info] = src.potential.Evaluate(rho, potreg, epsilon);
    relativeSteps = 10 .^ (-(2:8));
    errors = zeros(size(relativeSteps));
    localScale = max(abs(rho), max(info.sigma, 1e-8));
    for q = 1:numel(relativeSteps)
        step = relativeSteps(q) * localScale;
        plus = src.potential.Evaluate(rho + step, potreg, epsilon);
        minus = src.potential.Evaluate(rho - step, potreg, epsilon);
        finiteDifference = (plus - minus) ./ (2 * step);
        errors(q) = max(abs(finiteDifference - derivative));
    end
    bestErrors(j) = min(errors);
    assert(bestErrors(j) <= 1e-7, ...
        'Potential derivative error for %s is %.3e.', ...
        names{j}, bestErrors(j));
    if strcmp(names{j}, 'linear')
        assert(info.sigma == 0 && info.dp_at_zero == 1);
    else
        assert(info.sigma > 0 && info.dp_at_zero == 0);
    end
    fprintf('  potential derivative %-20s best max error %.3e\n', ...
        names{j}, bestErrors(j));
end
stats.names = names;
stats.best_max_error = bestErrors;
stats.maximum_error = max(bestErrors);
end
