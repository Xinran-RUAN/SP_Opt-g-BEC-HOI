% X is real-valued
function[E] = Obj_Func(X, data)
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;

%% get diff_X = DX 
D_X = fourier_diff(X, Lx);

%% kinetic energy
% % =================================================
% % 1: 1/8*int(|X'|^2/(X+vep)) 
% E_kin = h / 8 * sum(abs(D_X).^2 ./ (X + vep)); 
% % =================================================
% % 2: 1/8*int(|X'|^2/sqrt(X^2+vep^2)) 
E_kin = h / 8 * sum(abs(D_X).^2 ./ sqrt(X.^2 + vep^2)); 
% % =================================================
% % 3: 1/2*int(|Dx F|^2)
% [F_X, ~] = F_X_Func(X, vep);
% D_F_X = fourier_diff(F_X, Lambda);
% E_kin = h / 2 * sum(abs(D_F_X).^2);

%% potential energy
V = data.V;
E_pot = h * sum(V .* (sqrt(X.^2 + vep.^2) - vep));

%% beta - interaction
E_beta = 0.5 * beta * h * sum(X .^ 2);

%% delta - interaction
E_delta = 0.5 * delta * h * sum(abs(D_X).^2);

%% 高阶项 - 抑制高频部分
alpha = 1e-8;
Laplacian_X = fourier_diff_T(fourier_diff(X, Lx), Lx);
E_high = alpha * h * sum(Laplacian_X.^2);

%%
E = E_kin + E_pot + E_beta + E_delta + E_high;
