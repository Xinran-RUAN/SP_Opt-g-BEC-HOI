% 能量离散化 -- 傅立叶（DFT）拟谱离散
% 返回： E: 能量/目标函数
%       G: 梯度/能量变分, i.e.  G = delta E / delta X

function [E, G] = myDiscreteEnergy(X, data)
X = real(X); % check!
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;
N = data.Nx;
mu = 2 * pi / Lx * (-N/2:N/2-1)';
Dx = D_x(X, mu);
cDx = conjD_x(X, mu);
%% kinetic energy
E_kin = h / 8 * sum(abs(Dx).^2 ./ (X + vep)); 
G_kin = h * (tranD_x(cDx./(X+vep), mu) + conjtranD_x(Dx./(X+vep), mu) - abs(Dx).^2./(X+vep).^2);

%% potential energy
V = data.V;
E_pot = h * sum(V .* (sqrt(X.^2 + vep.^2) - vep));
G_pot = h * V .* X ./ sqrt(X.^2 + vep.^2);

%% beta - interaction
E_beta = 0.5 * beta * h * sum(X .^ 2);
G_beta = h * beta * X;

%% delta - interaction
E_delta = 0.5 * delta * h * sum(abs(Dx).^2);
G_delta = 0.5 * delta * h * (tranD_x(cDx, mu) + conjtranD_x(Dx, mu));

%%
E = E_kin + E_pot + E_beta + E_delta;
G = G_kin + G_pot + G_beta + G_delta;

%============================================================================
% kinetic energy term -- form 2
%   1/8*\frac{|\delta_x^+ \rho_j|^2}{\sqrt{\rho_{j+1/2}^2+\vep^2}}
%============================================================================
% % Boundary Condition (FD)
% X1=[0;X;0]; % X1(1)=X(1); X1(end)=X(end); % Neumann BC

% DX1 = (X1(2:end)-X1(1:end-1))/h; % \delta_x \rho_{j+1/2}
% AX1 = X1(1:end-1)/2+X1(2:end)/2;  % \rho_{j+1/2}
% X1_temp = DX1./(sqrt(AX1.^2+vep^2));
% X1_temp1 = X1_temp.^2.*AX1./(sqrt(AX1.^2+vep^2));
% 
% g1 = h/8*(-2*(X1_temp(2:end)-X1_temp(1:end-1))/h-0.5*(X1_temp1(1:end-1)+X1_temp1(2:end)));
% g = g1 + h*(V.*X./sqrt(X.^2+vep^2) + X*beta - delta*diff(X1,2)/h^2);
% f = h/8*sum(DX1.*X1_temp)+h*sum(V.*(sqrt(X.^2+vep^2)-vep)+X.*X*beta/2)+delta/2/h*(sum(diff(X1).^2));


% % % homo. Dirichlet BC (SP)
% lmd = pi/Lx*(1:N)';
% X_shift = 2/(N+1)*dst_iii(dst(X),N+1); % !!!: No ganruatee that X_shift>=0
% dX_shift = 2/(N+1)*dct_iii(lmd.*dst(X),N+1);
% F1_temp = dX_shift./(sqrt(X_shift.^2+vep^2));
% F2_temp = X_shift./(sqrt(X_shift.^2+vep^2));
% g1=h/4/(N+1)*(2*dst(lmd.*dct_ii(F1_temp,N+1))-dst(dst_ii(F2_temp.*F1_temp.^2,N+1))); % ... 
% % idst == 2/(N+1) * dst ??
% 
% G = g1+h*(V.*X./sqrt(X.^2+vep^2) + X*beta + delta*dst(lmd.^2.*idst(X)));
% E = E_kin + E_pot + E_beta + E_delta;


