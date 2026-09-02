clc;clear;
NN_ex=12;
beta=10; delta=10;  
mvep_ini=1e-1;
N_veplen = 5; NN_list_vep=[12,12,12,12,12,12,12];
E_test=zeros(N_veplen,1);

E_err=zeros(N_veplen,1);
l2_err=zeros(N_veplen,1);
h1_err=zeros(N_veplen,1);
linf_err=zeros(N_veplen,1);
% load exact solution
filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-1e-08-NN-',int2str(NN_ex),'_n.mat');
load(filename);N_length=length(Rho);
data_ex=data;
E_ex=E; % Rho_ex=Rho((N_length+3)/4:(3*N_length+1)/4); % x in [-32,32]
Rho_ex=Rho; 


% filename=strcat('GPE-SP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-NN-',int2str(NN_ex));
% load(filename)
% E_ex=f; Rho_ex=Phi.^2;

for mdvep=1:N_veplen
    NN=NN_list_vep(mdvep);    
    vep=mvep_ini/10^(mdvep-1);
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
%     disp(filename)
    load(filename); N_length=length(Rho);
%     if data.xmax == 32 && data_ex.xmax == 32
%         rho=finer(Rho((N_length+3)/4:(3*N_length+1)/4),NN_ex-NN); 
%     elseif data.xmax == 16 && data_ex.xmax == 32
%         rho=Rho;
%     else
%         disp('Error in the range of x')
%     end
    rho=Rho;
    
    rho_ex=Rho_ex; dh=1/2^NN_ex;
%     if NN>NN_ex
%         rho=Rho(1:2^(NN-NN_ex):end);
%         rho_ex=Rho_ex;
%         dh=1/2^NN_ex;
%     else
%         rho=Rho;
%         rho_ex=Rho_ex(1:2^(NN_ex-NN):end);
%         dh=1/2^NN;
%     end
    E_test(mdvep)=E;
    E_err(mdvep)=abs(E-E_ex);
    l2_err(mdvep)=sqrt(dh*sum((rho-rho_ex).^2));
%     h1_err(mdvep)=sqrt(sum(diff(rho-rho_ex).^2)/dh);
    linf_err(mdvep)=max(abs(rho-rho_ex));
    
    plot(rho-rho_ex)
end


E_err
l2_err
% h1_err
linf_err

% E_err(4:end) 
log(E_err(1:end-1)./E_err(2:end))/log(10)
log(l2_err(1:end-1)./l2_err(2:end))/log(10)
log(linf_err(1:end-1)./linf_err(2:end))/log(10)