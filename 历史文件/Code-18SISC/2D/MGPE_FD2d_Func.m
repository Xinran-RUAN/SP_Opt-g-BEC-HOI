function [f, g] = MGPE_FD2d_Func(data,X)

beta = data.beta;
delta = data.delta;
hx = data.dx; hy = data.dy;
V = data.V;
vep =data.vep;

X_MAT=reshape(X,data.Nx-2,data.Ny-2);
X_MAT_F = zeros(data.Nx,data.Ny);
X_MAT_F(2:end-1,2:end-1)=X_MAT;
DX1=(X_MAT_F(2:end,:)-X_MAT_F(1:end-1,:))/hx; % \delta_x \rho
DY1=(X_MAT_F(:,2:end)-X_MAT_F(:,1:end-1))/hy; % \delta_y \rho
D2X1=(X_MAT_F(1:end-2,:)-2*X_MAT_F(2:end-1,:)+X_MAT_F(3:end,:))/hx^2; % \delta_x^2 \rho
D2Y1=(X_MAT_F(:,1:end-2)-2*X_MAT_F(:,2:end-1)+X_MAT_F(:,3:end))/hy^2; % \delta_y^2 \rho
AX1=(abs(X_MAT_F(1:end-1,:))+abs(X_MAT_F(2:end,:)))/2; % \rho_{i+1/2,j}
AY1=(abs(X_MAT_F(:,1:end-1))+abs(X_MAT_F(:,2:end)))/2; % \rho_{i,j+1/2}
X1_temp=DX1./(AX1+vep); Y1_temp=DY1./(AY1+vep);


g1=hx*hy/8*(-2*(X1_temp(2:end,:)-X1_temp(1:end-1,:))/hx-0.5*(X1_temp(1:end-1,:).^2+X1_temp(2:end,:).^2).*sign(X_MAT_F(2:end-1,:)));
g2=hx*hy/8*(-2*(Y1_temp(:,2:end)-Y1_temp(:,1:end-1))/hy-0.5*(Y1_temp(:,1:end-1).^2+Y1_temp(:,2:end).^2).*sign(X_MAT_F(:,2:end-1)));
G = g1(:,2:end-1)+g2(2:end-1,:) + hx*hy*(V.*X_MAT./sqrt(X_MAT.^2+vep^2) + beta*X_MAT - delta*(D2X1(:,2:end-1)+D2Y1(2:end-1,:)));
g=reshape(G,(data.Nx-2)*(data.Ny-2),1);
f = hx*hy*sum(sum(V.*(sqrt(X_MAT.^2+vep^2)-vep)+beta/2*X_MAT.*X_MAT))+delta/2*hx*hy*(sum(sum(DX1.*DX1))+sum(sum(DY1.*DY1)))+hx*hy/8*(sum(sum(DX1.^2./(AX1+vep)))+sum(sum(DY1.^2./(AY1+vep))));

