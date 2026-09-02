% 输入：
%   初始 X；其他必要数据 data
% 程序中设置：
%   参数；
% 输出：
%   在可行域内寻找最小值点 Xc
% 同文件夹下需要如下函数：
%   1. 目标函数：Obj_Func
%   2. 目标函数梯度：GradObj_Func
%   3. 投影至某可行域的算法（计算p_L）：myP_L
function [Xc, k] = Routine_ISTA(X, data)
%===================================================
% 参数设置
%   L 由back-tracking的方式确定，初始值为L0，每次增大为eta*L0
%   当 \|X_c-X_o\| < tol 时算法停止，认为近似找到最优解
L0 = 1;
eta = 1.1;
tol = 1e-12; 
%===================================================
% 算法初始值设定
Xc = X;
k = 1;
err_l2 = 1;
max_iter = 1e3;
%===================================================
while(k==1 || err_l2 > tol)
    L = L0;
    %===================================================
    % 计算 Xn = p_L(Xc)
    grad_f = GradObj_Func(Xc, data);
    X_prox = pL_Func(Xc, L, data);
    % 使用back-tracking寻找L
    while(Obj_Func(X_prox, data) >  QL_Func(X_prox, Xc, L, data))  
        L = eta * L;
        X_prox = pL_Func(Xc, L, data);
    end
    % 更新
    Xo = Xc;
    Xc = X_prox; 
    k = k + 1;
    if k > max_iter
        break;
    end
    % L2-误差计算    
    err_l2 = sqrt(data.dx * sum((Xc - Xo).^2));
    disp(['ISTA:' ,num2str(err_l2)]);

    %% 测试能量下降速度
    % E(k) = Obj_Func(Xc, data);
    % plot(E(2:end)); title('Energy vs Iteration');
end
% hold off
