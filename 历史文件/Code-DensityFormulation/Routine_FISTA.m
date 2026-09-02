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
function [Xc, k] = Routine_FISTA(X, data)
%===================================================
% 参数设置
%   L初值为L0，每次以 eta*L0 的速度增大
%   当 \|X_c-X_o\| < tol 时算法停止，认为近似找到最优解
L0 = 5 * max(abs(GradObj_Func(X, data)));
eta = 1.1; 
tol = 1e-12; 
%===================================================
% 算法初始值设定
Xc = X;
Yc = Xc;
tc = 1;
k = 1;
err_l2 = 1;
max_iter = 1e5;
%===================================================
L = L0; % L似乎只需赋1个初值
while(k==1 || err_l2 > tol)    
    %===================================================
    Y_prox = pL_Func(Yc, L, data);
   
    % 使用back-tracking寻找L
    while(Obj_Func(Y_prox, data) > QL_Func(Y_prox, Yc, L, data))  
        L = eta * L;
        Y_prox = pL_Func(Yc, L, data);
    end
    % update
    Xo = Xc;
    Xc = Y_prox; 
    to = tc;
    tc = 0.5 * (1 + sqrt(1 + 4 * to^2));
    Yc = Xc + (to - 1) / tc * (Xc - Xo);
    k = k + 1;
    % 次数判断
    if k > max_iter
        break;
    end
    % L2-误差计算    
    err_l2 = sqrt(data.dx * sum((Xc - Xo).^2));
    disp(['FISTA:' ,num2str(err_l2), ', L = ', num2str(L)]);
end
