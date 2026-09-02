function expression = Expression(name)
%EXPRESSION Executable and readable definitions of supported p(rho).
%
% The function handles are the single numerical definition used by the
% potential energy, gradient, and Hessian curvature.  The text fields make
% the active model explicit in saved configurations and summaries.

name = lower(char(name));
switch name
    case 'linear'
        expression.value = @(rho, sigma) rho; %#ok<INUSD>
        expression.first_derivative = ...
            @(rho, sigma) ones(size(rho)); %#ok<INUSD>
        expression.second_derivative = ...
            @(rho, sigma) zeros(size(rho)); %#ok<INUSD>
        expression.text.sigma = '0';
        expression.text.value = 'rho';
        expression.text.first_derivative = '1';
        expression.text.second_derivative = '0';
    case {'sqrt_power', 'sqrt_same_scale', 'sqrt_squared_scale'}
        expression.value = @(rho, sigma) ...
            rho .^ 2 ./ (hypot(rho, sigma) + sigma);
        expression.first_derivative = @(rho, sigma) ...
            rho ./ hypot(rho, sigma);
        expression.second_derivative = @(rho, sigma) ...
            (sigma ./ hypot(rho, sigma)) .^ 2 ./ hypot(rho, sigma);
        expression.text.sigma = 'epsilon.^power';
        expression.text.value = ...
            'rho.^2 ./ (hypot(rho,sigma) + sigma)';
        expression.text.first_derivative = ...
            'rho ./ hypot(rho,sigma)';
        expression.text.second_derivative = ...
            '(sigma./hypot(rho,sigma)).^2 ./ hypot(rho,sigma)';
    otherwise
        error('src:potential:Expression:UnsupportedName', ...
            'Unsupported potential regularization "%s".', name);
end
expression.sigma_from_power = @(epsilon, power) epsilon .^ power;
end
