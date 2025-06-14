function [f, g] = MGPE_SP1d_Func(data,X)

beta = data.beta;
delta = data.delta;
h = data.dx;
V = data.V;
xmax=data.xmax;
xmin=data.xmin;
vep =data.vep;
N = data.Nx-2;
Lx=xmax-xmin;



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


% % homo. Dirichlet BC (SP)
lmd = pi/Lx*(1:N)';
X_shift = 2/(N+1)*dst_iii(dst(X),N+1); % !!!: No ganruatee that X_shift>=0
dX_shift = 2/(N+1)*dct_iii(lmd.*dst(X),N+1);
F1_temp = dX_shift./(sqrt(X_shift.^2+vep^2));
F2_temp = X_shift./(sqrt(X_shift.^2+vep^2));
g1=h/4/(N+1)*(2*dst(lmd.*dct_ii(F1_temp,N+1))-dst(dst_ii(F2_temp.*F1_temp.^2,N+1))); % ... 
% idst == 2/(N+1) * dst ??

g = g1+h*(V.*X./sqrt(X.^2+vep^2) + X*beta + delta*dst(lmd.^2.*idst(X)));
f = h/8*sum(dX_shift.*F1_temp)+h*sum(V.*(sqrt(X.^2+vep^2)-vep)+X.*X*beta/2+0.5*delta*X.*dst(lmd.^2.*idst(X)));


if nargout<3, return; end

% Chemical Potential
mu = f + h*sum(X.*X*beta/2+0.5*delta*X.*dst(lmd.^2.*idst(X)));