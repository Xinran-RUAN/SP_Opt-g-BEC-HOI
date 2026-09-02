function expression = KineticExpression(reg)
%KINETICEXPRESSION Executable/readable regularized Fisher expression.
%
% Derivative-based code obtains r, r', and r'' from EvaluateDenominator.
% This helper exposes the kinetic density itself as one direct expression,
% together with readable metadata for the selected denominator.

name = lower(char(reg.name));
expression.energy_density = @(Drho, denominator) ...
    Drho .^ 2 ./ (8 .* denominator);
expression.text.energy_density = '(Drho).^2 ./ (8 .* r_epsilon(rho))';
switch name
    case 'shift_smooth'
        expression.text.denominator = 'r_epsilon(rho) = rho + epsilon';
        expression.text.first_derivative = 'dr_epsilon(rho) = 1';
        expression.text.second_derivative = 'd2r_epsilon(rho) = 0';
    case {'piecewise_c1', 'piecewise_c2', 'piecewise_c3'}
        expression.text.denominator = ...
            ['r_epsilon(rho) = convexity-preserving ' name ' expression'];
        expression.text.first_derivative = ...
            'dr_epsilon(rho) from EvaluateDenominator';
        expression.text.second_derivative = ...
            'd2r_epsilon(rho) from EvaluateDenominator';
    otherwise
        error('src:regularization:KineticExpression:UnsupportedName', ...
            'Unsupported kinetic regularization "%s".', name);
end
end
