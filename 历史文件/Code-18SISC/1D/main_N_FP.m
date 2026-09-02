clc;clear;

beta=10; delta=100;  
mvep=1e-4; 
NN_ex=6; h_ex=1/NN_ex;
NN_list_N=[0,1,2,3,4];
NN_len_N = length(NN_list_N);
NN_end=NN_list_N(end);
E_test=zeros(NN_len_N,1);

E_err=zeros(NN_len_N,1);
l2_err=zeros(NN_len_N,1);
h1_err=zeros(NN_len_N,1);
linf_err=zeros(NN_len_N,1);
% load exact solution
filename=strcat('MGPE-FP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(mvep),'-NN-',int2str(NN_ex),'_n.mat');
load(filename)
E_ex=E; Rho_ex=Rho(1:end-1); % error in the boundary
x_ex=data.x;

for md_N=1:NN_len_N
    NN=NN_list_N(md_N);    
    filename=strcat('MGPE-FP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(mvep),'-NN-',int2str(NN),'_n.mat');
%     disp(filename)
    load(filename)
    x=data.x;
    Rho_test=finer_FP(Rho(1:end-1),NN_ex-NN); % Be careful: in whole space problem, boundary might be not 0

    E_err(md_N)=abs(E-E_ex);
    l2_err(md_N)=sqrt(sum(h_ex*(Rho_test-Rho_ex).^2));
    h1_err(md_N)=sqrt(sum(diff(Rho_test-Rho_ex).^2)/h_ex); % inaccurate for FP
    linf_err(md_N)=max(abs(Rho_test-Rho_ex));
        
    plot(Rho_test-Rho_ex); % hold on;
    pause(0.1);
end
% hold off;

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