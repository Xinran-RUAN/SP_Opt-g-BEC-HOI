function X0 = refine_3D(Rho)
[M,N,P]=size(Rho);
Rho_n=zeros(2*M-1,2*N-1,2*P-1);
Rho_n(1:2:end,1:2:end,1:2:end)=Rho;
Rho_n(2:2:end,1:2:end,1:2:end)=0.5*(Rho(1:end-1,:,:)+Rho(2:end,:,:));
Rho_n(1:2:end,2:2:end,1:2:end)=0.5*(Rho(:,1:end-1,:)+Rho(:,2:end,:));
Rho_n(1:2:end,1:2:end,2:2:end)=0.5*(Rho(:,:,1:end-1)+Rho(:,:,2:end));
Rho_n(2:2:end,2:2:end,1:2:end)=0.25*(Rho(1:end-1,1:end-1,:)+Rho(1:end-1,2:end,:)+Rho(2:end,1:end-1,:)+Rho(2:end,2:end,:));
Rho_n(1:2:end,2:2:end,2:2:end)=0.25*(Rho(:,1:end-1,1:end-1)+Rho(:,1:end-1,2:end)+Rho(:,2:end,1:end-1)+Rho(:,2:end,2:end));
Rho_n(2:2:end,1:2:end,2:2:end)=0.25*(Rho(1:end-1,:,1:end-1)+Rho(1:end-1,:,2:end)+Rho(2:end,:,1:end-1)+Rho(2:end,:,2:end));
Rho_n(2:2:end,2:2:end,2:2:end)=0.125*(Rho(1:end-1,1:end-1,1:end-1)+Rho(1:end-1,1:end-1,2:end)+Rho(2:end,1:end-1,1:end-1)+Rho(2:end,1:end-1,2:end)+...
    Rho(1:end-1,2:end,1:end-1)+Rho(1:end-1,2:end,2:end)+Rho(2:end,2:end,1:end-1)+Rho(2:end,2:end,2:end));
X0=Rho_n(2:end-1,2:end-1,2:end-1);