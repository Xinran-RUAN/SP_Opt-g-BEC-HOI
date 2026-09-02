% function MGPE_FP1d_Test1
% Change \gm if V is changed
% Qn: dh can' be too small
clear; clc; 
tic;
data.init = 'Gaussian'; 
Nbx = 5;
data.Potential = @(x) x.^2/2;
% data.Potential = @(x) x.^2/2;
data.beta = 1e1;
data.delta = 1e2;
data.xmin = -2^Nbx/2;
data.xmax = 2^Nbx/2;
NN_list = {0,1,2,3,4,5,6};
NN_start = NN_list{1}; 
NN = NN_start; 
NN_end = NN_list{end};
data.Nx = (data.xmax-data.xmin)*2^NN+1;

% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.x = (data.xmin:data.dx:data.xmax)';
% Potential function
x = data.x(1:end-1);
data.V = data.Potential(x);

vep_ini=1e-4;
veplen = 1;
Nlen = length(NN_list);

for dvep = 1:veplen 
    vep = vep_ini/2^(dvep-1);
    data.vep = vep;        
        for dN = 1:Nlen        
            NN = NN_list{dN}; data.Nx = (data.xmax-data.xmin)*2^NN+1;
            data.dx = (data.xmax-data.xmin)/(data.Nx-1);
            data.x = (data.xmin:data.dx:data.xmax)';
            % Potential function
            x = data.x(1:end-1); % start from 1 !!
            data.V = data.Potential(x);
            % Initial value
            if(dN == 1 && dvep == 1)  
%                 construct initial data
                X0 = MGPE_FP1d_Init(data); 
            elseif (dN == 1 && dvep > 1)
                Rho_n = Rho(1:2^(NN_end-dN+1):end);
                X0 = Rho_n(1:end-1);
            else
%                 X0 = zeros(2^(NN+Nbx),1);
%                 for Ni = 2:2:2^(NN+Nbx)-2
%                     X0(Ni) = 0.5*(Rho(Ni/2)+Rho(Ni/2+1));
%                 end
%                 for Ni = 1:2:2^(NN+Nbx)-1
%                     X0(Ni) = Rho((Ni+1)/2);
%                 end
                X0 = interpft(Rho(1:end-1),2^(NN+Nbx));
                % periodic BC
                
            end
            
            % update not only Rho, but also data by new values of NN
            [E,Rho] = MGPE_FP1d_Exam1(data,vep,NN,X0);           
            disp(E)
            
            filename = strcat('MGPE-FP1d-Bet-',int2str(data.beta),'-Del-',int2str(data.delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
            save(filename)
        end
      
end
toc
plot(data.x,Rho)
