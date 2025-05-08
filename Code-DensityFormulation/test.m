% 检验动能部分的梯度
% 当前计算方法：  dE_kin
% 原计算方法：    G_kin
% ============================================================
% 初始化data参数
vep = 1e-4; 
% 计算区域
Lx = 16;
h = 1 / 2^3;
data.xmax = Lx;
data.xmin = -Lx;
data.Nx = (data.xmax-data.xmin) / data.dx;
data.x = (data.xmin:data.dx:data.xmax-data.dx)';

% ============================================================
% 初始化数值解
x = data.x;
X = exp(-x.^2)/(pi^(1/2));

%% 
N = data.Nx;
Lambda = 2 * pi * 1i / Lx * [0:(N/2-1), (-N/2:-1)];
Lambda = reshape(Lambda, size(X)); % Lambda与X形状保持一致

%% 计算动能梯度
D_X = fourier_diff(X, Lambda);
G_X = D_X ./ (X + vep);
Dt_G_X = fourier_diff_T(G_X, Lambda);
dE_kin = h * (Dt_G_X / 4 - G_X.^2 / 8);

mu = 2 * pi / Lx * (-N/2:N/2-1)';
Dx = D_x(X, mu);
Dt_X = conjD_x(X, mu);
G_kin = h * (tranD_x(Dt_X./(X+vep), mu) + conjtranD_x(Dx./(X+vep), mu) - abs(Dx).^2./(X+vep).^2);

%% 通过对每个分量进行扰动，计算梯度的近似
Numerical_Grad = zeros(size(X));
nu = 1e-4;
for i = 1:length(X)
        e = zeros(size(X)); e(i) = 1;
        X_r = X + nu * e;
        X_l = X - nu * e;
        DX_r = fourier_diff(X_r, Lambda);
        DX_l = fourier_diff(X_l, Lambda);
        E_r = h / 8 * sum(abs(DX_r).^2 ./ (X_r + vep));
        E_l = h / 8 * sum(abs(DX_l).^2 ./ (X_l + vep));
        Numerical_Grad(i) = (E_r - E_l) / (2*nu);
end

%%
figure(1)
plot(Numerical_Grad, 'b-'); hold on;
plot(dE_kin, 'r--'); hold on;
plot(G_kin, 'k.'); hold off;

