% X is real-valued
function[dE] = GradObj_Func(X, data)
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;
V = data.V;

%% get diff_X = DX 
D_X = fourier_diff(X, Lx);

%% kinetic energy
% % =================================================
% % 1: 1/8*int(|X'|^2/(X+vep)) 
% G_X = D_X ./ (X + vep);
% Dt_G_X = fourier_diff_T(G_X, Lambda);
% dE_kin = h * (Dt_G_X / 4 - G_X.^2 / 8);
% % =================================================
% % 2: 1/8*int(|X'|^2/sqrt(X^2+vep^2)) 
G_X = D_X ./ sqrt(X.^2 + vep^2);
Dt_G_X = fourier_diff_T(G_X, Lx);
dE_kin = h * (Dt_G_X / 4 - G_X.^2 .* X ./ sqrt(X.^2 + vep^2) / 8);
% % =================================================
% [F_X, dF_dX] = F_X_Func(X, vep);
% D_F_X = fourier_diff(F_X, Lambda);
% Dt_D_F_X = fourier_diff_T(D_F_X, Lambda);
% dE_kin = h * dF_dX .* Dt_D_F_X;
 
%% potential energy
dE_pot = h * V .* X ./ sqrt(X.^2 + vep.^2);

%% beta - interaction
dE_beta = h * beta * X;

%% delta - interaction
Dt_D_X = fourier_diff_T(D_X, Lx);
dE_delta = h * delta * Dt_D_X;

%% 高阶项 - 抑制高频部分
alpha = 1e-8;
Laplacian_X = fourier_diff_T(fourier_diff(X, Lx), Lx);
biLaplacian_X = fourier_diff_T(fourier_diff(Laplacian_X, Lx), Lx);
dE_high = 2 * alpha * biLaplacian_X;

%%
dE = dE_kin + dE_pot + dE_beta + dE_delta + dE_high;
