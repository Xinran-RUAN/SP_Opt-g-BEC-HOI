function result = test_fisher_regularization_family()
%TEST_FISHER_REGULARIZATION_FAMILY Matching and convexity invariants.

epsilon = 1e-2;
ratios = [0; 1e-8; 0.1; 0.5; 0.999999; 1; 1.000001; 2; 10];
rho = epsilon * ratios;
families = cell(3, 1);
families{1}.r_epsilon = @(x) x + epsilon;
families{1}.dr_epsilon = @(x) ones(size(x));
families{1}.d2r_epsilon = @(x) zeros(size(x));
families{1}.epsilon = epsilon;
families{1}.name = 'shift';
families{1}.label = 'r_epsilon(rho)=rho+epsilon';
families{2} = src.regularization.PiecewiseFisherCm(epsilon, 1);
families{3} = src.regularization.PiecewiseFisherCm(epsilon, 2);

minimumR = inf;
minimumDr = inf;
maximumD2r = -inf;
for index = 1:numel(families)
    [families{index}, ~] = src.regularization.ValidateFisher( ...
        families{index}, rho);
    [r, dr, d2r] = src.regularization.EvaluateFisher( ...
        rho, families{index});
    minimumR = min(minimumR, min(r));
    minimumDr = min(minimumDr, min(dr));
    maximumD2r = max(maximumD2r, max(d2r));
    assert(all(r > 0), 'Fisher denominator must be positive.');
    assert(all(dr >= 0), 'Fisher denominator must be nondecreasing.');
    assert(all(d2r <= 1e-13), 'Fisher denominator must be concave.');
end

c1 = families{2};
c2 = families{3};
assert(abs(c1.r_epsilon(0) - epsilon / 2) <= 10 * eps(epsilon));
assert(abs(c2.r_epsilon(0) - 2 * epsilon / 3) <= 10 * eps(epsilon));

matchingStep = epsilon * 1e-13;
left = epsilon - matchingStep;
right = epsilon + matchingStep;
c1R = abs(c1.r_epsilon(left) - c1.r_epsilon(right)) ...
    / max(1, 2 * epsilon);
c1Dr = abs(c1.dr_epsilon(left) - c1.dr_epsilon(right));
c2R = abs(c2.r_epsilon(left) - c2.r_epsilon(right)) ...
    / max(1, 2 * epsilon);
c2Dr = abs(c2.dr_epsilon(left) - c2.dr_epsilon(right));
c2D2r = abs(c2.d2r_epsilon(left) - c2.d2r_epsilon(right)) ...
    / max(1, 1 / epsilon);
assert(c1R <= 1e-12 && c1Dr <= 1e-12, ...
    'Piecewise C1 matching failed.');
assert(c2R <= 1e-12 && c2Dr <= 1e-12 && c2D2r <= 1e-12, ...
    'Piecewise C2 matching failed.');

result.minimum_r = minimumR;
result.minimum_dr = minimumDr;
result.maximum_d2r = maximumD2r;
result.c1_matching_error = max(c1R, c1Dr);
result.c2_matching_error = max([c2R, c2Dr, c2D2r]);
fprintf(['test_fisher_regularization_family passed: C1 %.3e, ' ...
    'C2 %.3e, min(r) %.3e.\n'], result.c1_matching_error, ...
    result.c2_matching_error, result.minimum_r);
end
