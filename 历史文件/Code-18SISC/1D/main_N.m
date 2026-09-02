clc;clear;

beta=10; delta=10;  
mvep=1e-3; 
NN_ex=12; h_ex=1/NN_ex;
NN_list_N=[1,2,3,4,5,6];
NN_len_N = length(NN_list_N);
NN_end=NN_list_N(end);
E_test=zeros(NN_len_N,1);

E_err=zeros(NN_len_N,1);
l2_err=zeros(NN_len_N,1);
h1_err=zeros(NN_len_N,1);
linf_err=zeros(NN_len_N,1);
% load exact solution
% NN_ex=10;
% filename=strcat('GPE-SP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-NN-',int2str(NN_ex));
% load(filename)
% E_ex=f; Rho_ex=Phi.^2;
filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(mvep),'-NN-',int2str(NN_ex),'_n.mat');
load(filename)
E_ex=E; Rho_ex=Rho(2:end-1); % error in the boundary
x_ex=data.x;

for md_N=1:NN_len_N
    NN=NN_list_N(md_N);    
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(mvep),'-NN-',int2str(NN),'_n.mat');
%     disp(filename)
    load(filename)
    x=data.x;
%     Rho_test=finer(Rho,NN_ex-NN); % Be careful: in whole space problem, boundary might be not 0
    Rho_test=interp1(x(2:end-1),Rho(2:end-1),x_ex(2:end-1),'linear','extrap');
    E_err(md_N)=abs(E-E_ex);
    l2_err(md_N)=sqrt(sum(h_ex*(Rho_test-Rho_ex).^2));
    h1_err(md_N)=sqrt(sum(diff(Rho_test-Rho_ex).^2)/h_ex);
    linf_err(md_N)=max(abs(Rho_test-Rho_ex));
    
%     if NN>NN_ex
%         rho=Rho(1:2^(NN-NN_ex):end);
%         rho_ex=Rho_ex;
%         dh=1/2^NN_ex;
%     else
%         rho=Rho;
%         rho_ex=Rho_ex(1:2^(NN_ex-NN):end);
%         dh=1/2^NN;
%     end
%     E_test(md_N)=E;
%     E_err(md_N)=abs(E-E_ex);
%     l2_err(md_N)=sqrt(dh*sum((rho-rho_ex).^2));
%     h1_err(md_N)=sqrt(sum((diff(rho)-diff(rho_ex)).^2)/dh);
% %     h1_err(md_N)=sqrt(sum(diff(rho-rho_ex).^2)/dh);
%     linf_err(md_N)=max(abs(rho-rho_ex));
%     
    plot(Rho_test-Rho_ex)
end

format long
E_err
l2_err
h1_err
linf_err

% E_err(4:end) 
log(E_err(1:end-1)./E_err(2:end))/log(2)
log(l2_err(1:end-1)./l2_err(2:end))/log(2)
log(h1_err(1:end-1)./h1_err(2:end))/log(2)
log(linf_err(1:end-1)./linf_err(2:end))/log(2)