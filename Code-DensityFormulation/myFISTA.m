% FISTA -- 未检查
function [Xc] = myFISTA(X, data)
% parameters
L0 = 1;
eta = 1.1;
tol = 1e-12;
h = data.dx;
% initial data
Xc = X;
y = Xc;
k = 1;
t = 1;
err = 1;

%
while(k==1 || err > tol)
    G = GradObj_Func(Xc, data);
    L = ceil(max(abs(G)));
    X_temp = myP_L(y, L, data);
    E = myDiscreteEnergy(X_temp, data);
    %
    QL_1 = myQ_L(X_temp, y, L, data);
    % % for test
    QL = QL_Func(X_temp, y, L, @Obj_Func, @GradObj_Func, data);
    disp(abs(QL_1-QL));

    while(real(E) > real(QL))   % check!
        L = eta * L;
        X_temp = myP_L(y, L, data);
        if max(isnan(X_temp)) == 1
            continue;
        end
        E = myDiscreteEnergy(X_temp, data);
        QL = myQ_L(X_temp, y, L, data);
    end
    % update
    Xo = Xc;
    Xc = X_temp; 
    t_temp = 0.5 * (1 + sqrt(1 + 4 * t^2));
    y = Xc + (t - 1) / t_temp * (Xc - Xo);
    t = t_temp;
    err = sqrt(h * sum(abs(Xc-Xo).^2));
    k = k + 1;
    figure(1)
    plot(data.x, real(Xc));
    title(strcat('k=',num2str(k), ', error=',num2str(err), ', L=',num2str(L)));

end
