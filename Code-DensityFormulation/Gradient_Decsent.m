% 梯度下降法
% 通过控制步长，保证 min(Xn)的最小值不可离0太远
function [Xn] = Gradient_Decsent(X, data)
[~, G] = myDiscreteEnergy(X, data);
tol = -1e-3;
min_X = -1;
t = 1;
while min_X < tol
    Xn = X - t * G; 
    min_X = min(Xn);
    t = t / 2;
end
