 % function MGPE_FD1d_Test1
 % do this if not have been done: addpath /Users/ruanxinran/Desktop/matlab_bin/TFOCS-master
% ATTENTION: in the file MGPE_FD1d_Init.m, \gm is used if init='TFA_D'
% Change \gm if V is changed
% Qn: dh can' be too small
clear; clc; 
tic;
% data.init = 'TFA_D'; 
data.init = 'Gaussian'; 
Nbx=5;
% Nbx=1;
data.Potential = @(x) x.^2/2;
data.beta =10;
data.delta = 1e4;
data.xmin = -2^Nbx/2;
data.xmax = 2^Nbx/2;
NN_list={0,1,2,3,4,5,6};
% NN_list={6};
NN_start=NN_list{1}; NN=NN_start; NN_end=NN_list{end};
data.Nx = (data.xmax-data.xmin)*2^NN+1;

% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.x = (data.xmin:data.dx:data.xmax)';
% Potential function
x = data.x(2:end-1);
data.V = data.Potential(x);

vep_ini=1e-8;
veplen = 1; 
Nlen = length(NN_list);



for dvep=1:veplen 
    vep=vep_ini/10^(dvep-1);
    data.vep=vep;
    NN = NN_list{1};
    data.Nx = (data.xmax-data.xmin)*2^NN+1;
    data.dx = (data.xmax-data.xmin)/(data.Nx-1);
    data.x = (data.xmin:data.dx:data.xmax)';
    % Potential function
    x = data.x(2:end-1);
    data.V = data.Potential(x);
   
    % Initial value    
    if(dvep==1)  
        % construct initial data
        X0=MGPE_FD1d_Init(data); 
    else 
        X0=Rho(2:end-1);
    end
        
    [E,Rho] = MGPE_FD1d_Exam1(data,vep,NN,X0);
    
end

for dN = 1:Nlen        
    NN = NN_list{dN}; data.Nx = (data.xmax-data.xmin)*2^NN+1;
    data.dx = (data.xmax-data.xmin)/(data.Nx-1);
    data.x = (data.xmin:data.dx:data.xmax)';
    % Potential function
    x = data.x(2:end-1);
    data.V = data.Potential(x);
    % Initial value
    if(dN==1)  
        X0=Rho(2:end-1);
    else
        X0=zeros(2^(NN+Nbx)-1,1);
        for Ni=3:2:2^(NN+Nbx)-3
            X0(Ni)=0.5*(Rho((Ni+1)/2+1)+Rho((Ni-1)/2+1));
        end
        for Ni=2:2:2^(NN+Nbx)-2
            X0(Ni)=Rho(Ni/2+1);
        end
        % zero BC
                X0(1)=0.5*Rho(2);
                X0(end)=0.5*Rho(end-1);
        % Neumann BC
%         X0(1)=Rho(2);
%         X0(end)=Rho(end-1);
    end

    % update not only Rho, but also data by new values of NN
    [E,Rho] = MGPE_FD1d_Exam1(data,vep,NN,X0);

    disp(E)
%             plot(data.x,Rho); hold on;

    filename=strcat('MGPE-FD1d-Bet-',int2str(data.beta),'-Del-',int2str(data.delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
%             filename=strcat('MGPE-FD1d-Bet-',int2str(data.beta),'-Del-',int2str(data.delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_BOX_n.mat');
    save(filename)
end
toc
plot(data.x,Rho)
