% 检验动能部分的梯度
% 当前计算方法：  dE_kin
% 原计算方法：    G_kin
% ============================================================
% 初始化data参数
vep = 1e-8; 
% 计算区域
Lx = 16;
h = 1 / 2^3;
data.dx = h;
data.xmax = Lx;
data.xmin = -Lx;
data.Nx = (data.xmax-data.xmin) / data.dx;
data.x = (data.xmin:data.dx:data.xmax-data.dx)';

% ============================================================
% 初始化数值解
x = data.x;
X = exp(-x.^2)/(pi^(1/2));
X_in = X(2:end-1);

%% 
N = data.Nx;
Lambda = 2 * pi * 1i / (2* Lx) * [0:(N/2-1), (-N/2:-1)];
Lambda = reshape(Lambda, size(X)); % Lambda与X形状保持一致

%% 计算动能梯度
% =====================================================
% D_X = fourier_diff(X, Lambda);
% G_X = D_X ./ (X + vep);
% Dt_G_X = fourier_diff_T(G_X, Lambda);
% dE_kin = h * (Dt_G_X / 4 - G_X.^2 / 8);
% 
% mu = 2 * pi / Lx * (-N/2:N/2-1)';
% Dx = D_x(X, mu);
% Dt_X = conjD_x(X, mu);
% G_kin = h * (tranD_x(Dt_X./(X+vep), mu) + conjtranD_x(Dx./(X+vep), mu) - abs(Dx).^2./(X+vep).^2);
% =====================================================
% [F_X, dF_dX] = F_X_Func(X, vep);
% D_F_X = fourier_diff(F_X, Lambda);
% Dt_D_F_X = fourier_diff_T(D_F_X, Lambda);
% dE_kin = h * dF_dX .* Dt_D_F_X;
% =====================================================
% [F_X, dF_dX] = F_X_Func(X_in, vep);
% D_F_X = DST_diff(F_X, 2 * Lx);
% Dt_D_F_X = DST_diff_T(D_F_X, 2 * Lx);
% dE_kin = h * dF_dX .* Dt_D_F_X;

%% 通过对每个分量进行扰动，计算梯度的近似
Numerical_Grad = zeros(size(X));
nu = 1e-14;
for i = 1:length(X)
        e = zeros(size(X)); e(i) = 1;
        X_r = X + nu * e;
        X_l = X - nu * e;
        % =====================================================
        % DX_r = fourier_diff(X_r, Lambda);
        % DX_l = fourier_diff(X_l, Lambda);
        % E_r = h / 8 * sum(abs(DX_r).^2 ./ (X_r + vep));
        % E_l = h / 8 * sum(abs(DX_l).^2 ./ (X_l + vep));
        % =====================================================
        [F_X_r, ~] = F_X_Func(X_r, vep);
        D_F_X_r = fourier_diff(F_X_r, Lambda);
        E_r_f = h / 2 * sum(abs(D_F_X_r).^2);
        [F_X_l, ~] = F_X_Func(X_l, vep);
        D_F_X_l = fourier_diff(F_X_l, Lambda);
        E_l_f = h / 2 * sum(abs(D_F_X_l).^2);
        % =====================================================
        Numerical_Grad(i) = (E_r - E_l) / (2*nu);
end

for i = 1:length(X_in)
        % =====================================================
        e = zeros(size(X_in)); e(i) = 1;
        X_in_r = X_in + nu * e;
        X_in_l = X_in - nu * e;
        [F_X_in_r, ~] = F_X_Func(X_in_r, vep);
        D_F_X_in_r = DST_diff(F_X_in_r, 2 * Lx);
        E_r = h / 2 * sum(abs(D_F_X_in_r(2:end-1)).^2)+h/4*abs(D_F_X_in_r(1)).^2+h/4*abs(D_F_X_in_r(end)).^2;
        [F_X_in_l, ~] = F_X_Func(X_in_l, vep);
        D_F_X_in_l = DST_diff(F_X_in_l, 2 * Lx);
        E_l = h / 2 * sum(abs(D_F_X_in_l).^2)+h/4*abs(D_F_X_in_l(1)).^2+h/4*abs(D_F_X_in_l(end)).^2;
        % =====================================================
        Numerical_Grad(i) = (E_r - E_l) / (2*nu);
end

%%
% figure(1)
% plot(Numerical_Grad, 'b-'); hold on;
% plot(dE_kin, 'r--'); hold off;
% % plot(G_kin, 'k.'); hold off;

% compare derivatives
figure(2)
plot(D_F_X_in_l); hold on;
plot(-2*x .* exp(-x.^2)/(pi^(1/2)), 'k'); hold on; % 根据具体例子修改表达式
plot(D_F_X_l, 'r--'); hold off;
