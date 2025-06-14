function [f, g, mu] = MGPE_FD3d_Func( data,X)

beta = data.beta;
delta = data.delta;
hx = data.dx; hy = data.dy; hz = data.dz;
V = data.V;
vep =data.vep;

X_MAT=reshape(X,data.Nx-2,data.Ny-2,data.Nz-2);
X_MAT_F = zeros(data.Nx,data.Ny,data.Nz);
X_MAT_F(2:end-1,2:end-1,2:end-1)=X_MAT;
DX1=(X_MAT_F(2:end,:,:)-X_MAT_F(1:end-1,:,:))/hx; % \delta_x \rho
DY1=(X_MAT_F(:,2:end,:)-X_MAT_F(:,1:end-1,:))/hy; % \delta_y \rho
DZ1=(X_MAT_F(:,:,2:end)-X_MAT_F(:,:,1:end-1))/hz; % \delta_z \rho
D2X1=(X_MAT_F(1:end-2,:,:)-2*X_MAT_F(2:end-1,:,:)+X_MAT_F(3:end,:,:))/hx^2; % \delta_x^2 \rho
D2Y1=(X_MAT_F(:,1:end-2,:)-2*X_MAT_F(:,2:end-1,:)+X_MAT_F(:,3:end,:))/hy^2; % \delta_y^2 \rho
D2Z1=(X_MAT_F(:,:,1:end-2)-2*X_MAT_F(:,:,2:end-1)+X_MAT_F(:,:,3:end))/hz^2; % \delta_z^2 \rho
AX1=(abs(X_MAT_F(1:end-1,:,:))+abs(X_MAT_F(2:end,:,:)))/2; % \rho_{i+1/2,j,k}
AY1=(abs(X_MAT_F(:,1:end-1,:))+abs(X_MAT_F(:,2:end,:)))/2; % \rho_{i,j+1/2,k}
AZ1=(abs(X_MAT_F(:,:,1:end-1))+abs(X_MAT_F(:,:,2:end)))/2; % \rho_{i,j,k+1/2}
X1_temp=DX1./(AX1+vep); Y1_temp=DY1./(AY1+vep); Z1_temp=DZ1./(AZ1+vep);


g1=hx*hy*hz/8*(-2*(X1_temp(2:end,:,:)-X1_temp(1:end-1,:,:))/hx-0.5*(X1_temp(1:end-1,:,:).^2+X1_temp(2:end,:,:).^2).*sign(X_MAT_F(2:end-1,:,:)));
g2=hx*hy*hz/8*(-2*(Y1_temp(:,2:end,:)-Y1_temp(:,1:end-1,:))/hy-0.5*(Y1_temp(:,1:end-1,:).^2+Y1_temp(:,2:end,:).^2).*sign(X_MAT_F(:,2:end-1,:)));
g3=hx*hy*hz/8*(-2*(Z1_temp(:,:,2:end)-Z1_temp(:,:,1:end-1))/hz-0.5*(Z1_temp(:,:,1:end-1).^2+Z1_temp(:,:,2:end).^2).*sign(X_MAT_F(:,:,2:end-1)));
G = g1(:,2:end-1,2:end-1)+g2(2:end-1,:,2:end-1)+g3(2:end-1,2:end-1,:) + hx*hy*hz*(V.*X_MAT./sqrt(X_MAT.^2+vep^2) + beta*X_MAT - delta*(D2X1(:,2:end-1,2:end-1)+D2Y1(2:end-1,:,2:end-1)+D2Z1(2:end-1,2:end-1,:)));
g=reshape(G,(data.Nx-2)*(data.Ny-2)*(data.Nz-2),1);

f = hx*hy*hz*sum(sum(sum(V.*(sqrt(X_MAT.^2+vep^2)-vep)+beta/2*X_MAT.*X_MAT)))+delta/2*hx*hy*hz*(sum(sum(sum(DX1.*DX1)))+sum(sum(sum(DY1.*DY1)))+sum(sum(sum(DZ1.*DZ1))))...
    +hx*hy*hz/8*(sum(sum(sum(DX1.^2./(AX1+vep))))+sum(sum(sum(DY1.^2./(AY1+vep))))+sum(sum(sum(DZ1.^2./(AZ1+vep)))));
