 % function MGPE_FD2d_Test1
%  addpath /Users/ruanxinran/Desktop/matlab_bin/TFOCS-master !!!
% ATTENTION: in the file MGPE_FD1d_Init.m, \gm is used if init='TFA_D'
% Change \gm if V is changed
% Qn: dh can' be too small
clear; tic; 
data.init = 'Gaussian'; 
Nbx=1; Nby=Nbx;
data.Potential = @(x,y) 0*(x.^2+y.^2)/2;
data.beta = 1;
data.delta = 1;
data.xmin = -2^Nbx/2; data.xmax = 2^Nbx/2;
data.ymin = -2^Nby/2; data.ymax = 2^Nby/2;
NN_list={1,2,3,4,5,6};
NN_start=NN_list{1}; NN=NN_start; NN_end=NN_list{end};
data.Nx = (data.xmax-data.xmin)*2^NN+1;
data.Ny = (data.ymax-data.ymin)*2^NN+1;

% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.dy = (data.ymax-data.ymin)/(data.Ny-1);
data.x = (data.xmin:data.dx:data.xmax)';
data.y = (data.ymin:data.dy:data.ymax)';
% Potential function
x = data.x(2:end-1); y = data.y(2:end-1);
[X,Y]=ndgrid(x,y);
data.V = data.Potential(X,Y);

vep_ini=1e-4;
veplen = 1;
Nlen = length(NN_list);

 



for dvep=1:veplen 
    vep=vep_ini/2^(dvep-1);
    data.vep=vep;
   
        
        for dN = 1:Nlen        
            NN = NN_list{dN}; 
            data.Nx = (data.xmax-data.xmin)*2^NN+1;
            data.Ny = (data.ymax-data.ymin)*2^NN+1;
            data.dx = (data.xmax-data.xmin)/(data.Nx-1);
            data.dy = (data.ymax-data.ymin)/(data.Ny-1);
            data.x = (data.xmin:data.dx:data.xmax)';
            data.y = (data.ymin:data.dy:data.ymax)';
            % Potential function
            x = data.x(2:end-1); y = data.y(2:end-1);
            [X,Y]=ndgrid(x,y);
            data.V = data.Potential(X,Y);
            % Initial value
            if(dN==1 && dvep==1)   
                X0=MGPE_FD2d_Init(data); 
            elseif (dN==1 && dvep>1)
                Rho_n=Rho(1:2^(NN_end-dN+1):end,1:2^(NN_end-dN+1):end);
                X0=Rho_n(2:end-1,2:end-1);
            else
                X0=refine_2D(Rho);            
            end
            % update not only Rho, but also data by new values of NN
            [E,Rho] = MGPE_FD2d_Exam1(data,vep,NN,X0);
%             plot(data.x,Rho); hold on;
            % -------test-------------------------------------------------
    %             Rho_t=fmincon(fun,X0,A,b,Aeq,beq,lb,ub); Rho_t1=[0;Rho_t;0];
    %             plot(data.x,Rho_t1,'r--'); hold off;
    %             f1=fun(Rho_t);
            %-------------------------------------------------------------

%             filename=strcat('MGPE-FD1d-Bet-',int2str(data.beta),'-Del-',int2str(data.delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
%             save(filename)
        end
    
%         X0=Rho(2:end-1);
%         [E,mu,Rho,data] = MGPE_FD1d_Exam1(data,vep,NN,X0,solver,opts);
%         filename=strcat('MGPE-FD1d-Bet-',int2str(data.beta),'-Del-',int2str(data.delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
%         save(filename)
    
    


    
end
toc


%% save .mat
fname = strcat('MGPE_OptRegE_2D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep),'_Box'); % Har: harmonic potential Opt: optical potential
save(strcat(fname,'.mat'),'E','Rho','data','vep');


%% 2-D plot
textsize=30;
vep=1e-4; % change by hand in fname
data.beta=100; data.delta=500; 
fname = strcat('MGPE_OptRegE_2D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep),'_Box');
fname_fig = strcat('MGPE_OptRegE_2D_',num2str(data.beta),'_',num2str(data.delta),'_1e-4_Box');
filename=strcat(fname,'.mat');
load(filename)
fig=figure;
image(data.x,data.y,Rho,'CDataMapping','scaled')
colormap(jet)
title(strcat('(d) $\beta=',int2str(data.beta),' ,\delta=',int2str(data.delta),'$'),'Interpreter','latex');
set(gca,'xtick',[data.xmin,0,data.xmax],'ytick',[data.ymin,0,data.ymax],'YDir','normal');
set(gca,'FontSize',textsize);
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
% xlim([-2,2]); ylim([-2,2]);
axis equal; axis tight; 
set(gca,'YTick',[-1,0,1],'XTick',[-1,0,1])
caxis([0,0.7]);
cbh=colorbar; set(cbh,'YTick',[0:0.2:0.7]); % colorbar
saveas(fig,strcat(fname_fig,'.fig'));
print(fig,'-depsc',strcat(fname_fig,'.eps'));