function [f, g] = MGPE_FD1d_Func(data,X)

beta = data.beta;
delta = data.delta;
h = data.dx;
V = data.V;
vep =data.vep;


X1=[0;X;0]; % X1(1)=X(1); X1(end)=X(end); % Neumann BC


%============================================================================
% kinetic energy term -- form 1
%   1/4*\frac{|\delta_x^+ \rho_j|^2}{|\rho_j|+|\rho_{j+1}|+2\vep}
%============================================================================
% DX1 = (X1(2:end)-X1(1:end-1))/h; % \delta_x \rho_{j+1/2}
% AX1 = abs(X1(1:end-1))/2+abs(X1(2:end))/2+vep;  % \rho_{j+1/2}+\vep
% X1_temp = DX1./(AX1);
% 
% g1 = h/8*(-2*(X1_temp(2:end)-X1_temp(1:end-1))/h-0.5*(X1_temp(1:end-1).^2+X1_temp(2:end).^2).*sign(X));
% g = g1 + h*(V.*X./sqrt(X.^2+vep^2) + X*beta - delta*diff(X1,2)/h^2);
% f = h*sum(V.*(sqrt(X.^2+vep^2)-vep)+X.*X*beta/2)+delta/2/h*(sum(diff(X1).^2))+h/8*sum(DX1.^2./(AX1));


%============================================================================
% kinetic energy term -- form 2
%   1/8*\frac{|\delta_x^+ \rho_j|^2}{\sqrt{\rho_{j+1/2}^2+\vep^2}}
%============================================================================
DX1 = (X1(2:end)-X1(1:end-1))/h; % \delta_x \rho_{j+1/2}
AX1 = X1(1:end-1)/2+X1(2:end)/2;  % \rho_{j+1/2}
X1_temp = DX1./(sqrt(AX1.^2+vep^2));
X1_temp1 = X1_temp.^2.*AX1./(sqrt(AX1.^2+vep^2));

g1 = h/8*(-2*(X1_temp(2:end)-X1_temp(1:end-1))/h-0.5*(X1_temp1(1:end-1)+X1_temp1(2:end)));
g = g1 + h*(V.*X./sqrt(X.^2+vep^2) + X*beta - delta*diff(X1,2)/h^2);
f = h*sum(V.*(sqrt(X.^2+vep^2)-vep)+X.*X*beta/2)+delta/2/h*(sum(diff(X1).^2))+h/8*sum(DX1.*X1_temp);
% Attention: the energy computed in form 2 is larger than the energy
% computed in form 1 since AX1+vep>sqrt{AX1^2+vep^2}

% A potential problem: less regular than form 1??
% dut to the cut-off error when vep<<1 and rho not near 0??
