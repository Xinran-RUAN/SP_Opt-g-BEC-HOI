% Compare 
clc; clear;
beta=10; delta=10; NN_ex=8; vep_tt=1e-6;
filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep_tt),'-NN-',int2str(NN_ex),'_n.mat');
load(filename); x_ex_full=data.x;
x_ex=x; h_ex=x(2)-x(1); Rho_ex=Rho; N_Rho=length(Rho); E_ex=E; 
NN_test={1,2,3,4,5,6,7}; N=length(NN_test);
RHO=zeros(N_Rho,N);
l2_err=zeros(N,1);h1_err=zeros(N,1);
E_err=zeros(N,1);
linf_err=zeros(N,1);
Rho_test=zeros(N_Rho,1);

for jj=1:N
    N_jj=NN_test{jj};
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep_tt),'-NN-',int2str(N_jj),'_n.mat');
    load(filename);
    Rho_test=finer(Rho,NN_ex-N_jj);
    RHO(:,jj)=Rho_test; h=x(2)-x(1);
    E_err(jj)=abs(E-E_ex);
    l2_err(jj)=sqrt(sum(h_ex*(Rho_test-Rho_ex).^2));
%     h1_err(jj)=sqrt(sum(h_ex*diff(Rho_test-Rho_ex).^2));
    h1_err(jj)=sqrt(sum(diff(Rho_test-Rho_ex).^2)/h_ex);
    linf_err(jj)=max(abs(Rho_test-Rho_ex));
%     plot(x_ex_full,Rho_test-Rho_ex); hold on;
end
l2_err
% log(l2_err(1:end-1)./l2_err(2:end))/log(2)
h1_err
% log(h1_err(1:end-1)./h1_err(2:end))/log(2)
linf_err
% log(linf_err(1:end-1)./linf_err(2:end))/log(2)
E_err
% log(E_err(1:end-1)./E_err(2:end))/log(2)
% hold off;