function [g] = MGPE_FD2d_Grad(X, opt)
data=opt.data;
%=======================================================
% Gradient for MGPE
%=======================================================
beta = data.beta;
delta = data.delta;
hx = data.dx; hy = data.dy;
V = data.V;
vep =data.vep;

X_MAT=reshape(X,data.Nx-2,data.Ny-2);
X_MAT_F = zeros(data.Nx,data.Ny); X_MAT_F(2:end-1,2:end-1)=X_MAT;
DX1=(X_MAT_F(2:end,:)-X_MAT_F(1:end-1,:))/hx; % \delta_x \rho
DY1=(X_MAT_F(:,2:end)-X_MAT_F(:,1:end-1))/hy; % \delta_y \rho
D2X1=(X_MAT_F(1:end-2,:)-2*X_MAT_F(2:end-1,:)+X_MAT_F(3:end,:))/hx^2; % \delta_x^2 \rho
D2Y1=(X_MAT_F(:,1:end-2)-2*X_MAT_F(:,2:end-1)+X_MAT_F(:,3:end))/hy^2; % \delta_y^2 \rho
AX1=(X_MAT_F(1:end-1,:)+X_MAT_F(2:end,:))/2; % \rho_{i+1/2,j}
AY1=(X_MAT_F(:,1:end-1)+X_MAT_F(:,2:end))/2; % \rho_{i,j+1/2}
X1_temp=DX1./(AX1+vep); Y1_temp=DY1./(AY1+vep);


g1=hx*hy/8*(-2*(X1_temp(2:end,:)-X1_temp(1:end-1,:))/hx-0.5*(X1_temp(1:end-1,:).^2+X1_temp(2:end,:).^2));
g2=hx*hy/8*(-2*(Y1_temp(:,2:end)-Y1_temp(:,1:end-1))/hy-0.5*(Y1_temp(:,1:end-1).^2+Y1_temp(:,2:end).^2));
G = g1(:,2:end-1)+g2(2:end-1,:) + hx*hy*(V + beta*X_MAT - delta*(D2X1(:,2:end-1)+D2Y1(2:end-1,:)));

g=reshape(G,(data.Nx-2)*(data.Ny-2),1);