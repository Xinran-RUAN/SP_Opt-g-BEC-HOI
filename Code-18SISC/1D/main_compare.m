clc;clear;
beta_list=[10,20,40,80,1000]; N_beta=length(beta_list);
delta=10; 
NN=8; Ndx=1/2^NN;
Linewidth=2;textsize=25;
vep=1e-4;
% Linestyle={'-bs','-rd','-mo'};
Linestyle={'--','-.','.'};
l2_err=zeros(N_beta,1);

for jj=1:N_beta
    beta=beta_list(jj);
    % load exact solution
    filename=strcat('GPE-SP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-NN-6');
    load(filename)
    x_ex=data.x; dx_ex=x_ex(2)-x_ex(1);
    Rho_ex=Phi.^2; 

    % load rho_g^\vep
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
    load(filename) 
    x=data.x;
%     Rho_t=Rho(1:2^(NN-6):end);
    Rho_t=interp1(x,Rho,x_ex);
    
    l2_err(jj)=sqrt(dx_ex*sum((Rho_t-Rho_ex).^2));
    
    plot(Rho_t-Rho_ex)
    hold on;
    
end
hold off;
l2_err

