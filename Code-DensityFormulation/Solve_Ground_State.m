% 
function [E, Xn] = Solve_Ground_State(X0, data)
% 
err = 1;
tol = 1e-12;
Xn = X0;
while err > tol
    Xo = Xn;
    %% 梯度下降法
    % mass = 1;
    % domain_size = data.xmax - data.xmin;
    % Xn = Gradient_Decsent(Xo, data);  % 梯度下降
    % Xn = myPosConProj(Xn, mass, domain_size, 1e-12); % 投影至可行集： mass = dx * sum(x), dx = domain_size / Nx
    %% FISTA
    % Xn = myFISTA(Xo, data);
    Xn = Routine_FISTA(Xo, data);
    %%
    % err = max(abs(Xn - Xo));
    err = sqrt(data.dx * sum((Xn - Xo).^2));
    disp(err);
end

% Print Energy and Chemical Potential
[E, ~] = myDiscreteEnergy(Xn, data);
