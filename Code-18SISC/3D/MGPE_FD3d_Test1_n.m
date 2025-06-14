 % function MGPE_FD3d_Test1
 % addpath /Users/ruanxinran/Desktop/matlab_bin/TFOCS-master !!!
% ATTENTION: in the file MGPE_FD1d_Init.m, \gm is used if init='TFA_D'
% Change \gm if V is changed
% Qn: dh can' be too small
clear; tic;
data.init = 'Gaussian'; Pot='Opt';
Nbx=5; Nby=3; Nbz=3;
% data.Potential = @(x,y,z) (x.^2+4*y.^2+4*z.^2)/2+(2*sin(x).^2+10*sin(2*y).^2+20*sin(3*z).^2); % Opt1
data.Potential = @(x,y,z) (x.^2+4*y.^2+4*z.^2)/2+20*(sin(x).^2+sin(y).^2+sin(z).^2); % Opt
% data.Potential = @(x,y,z) (x.^2+4*y.^2+4*z.^2)/2+100*(sin(x).^2+sin(y).^2+sin(z).^2); % Opt2
% data.Potential = @(x,y,z) (x.^2+4*y.^2+4*z.^2)/2+40*(sin(x).^2+sin(y).^2+sin(z).^2); % Opt3
% data.Potential = @(x,y,z) (x.^2+4*y.^2+4*z.^2)/2; % Har
% data.Potential = @(x,y,z) 0*(x.^2+4*y.^2+4*z.^2)/2; % Box
data.beta = 1;
data.delta = 1;
data.xmin = -2^Nbx/2; data.xmax = 2^Nbx/2;
data.ymin = -2^Nby/2; data.ymax = 2^Nby/2;
data.zmin = -2^Nbz/2; data.zmax = 2^Nbz/2;
% NN_list={1,2,3,4,5};% box
NN_list={1,2,3};% har opt
NN_start=NN_list{1}; NN=NN_start; NN_end=NN_list{end};
data.Nx = (data.xmax-data.xmin)*2^NN+1;
data.Ny = (data.ymax-data.ymin)*2^NN+1;
data.Nz = (data.zmax-data.zmin)*2^NN+1;

% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.dy = (data.ymax-data.ymin)/(data.Ny-1);
data.dz = (data.zmax-data.zmin)/(data.Nz-1);
data.x = (data.xmin:data.dx:data.xmax)';
data.y = (data.ymin:data.dy:data.ymax)';
data.z = (data.zmin:data.dz:data.zmax)';
% Potential function
x = data.x(2:end-1); y = data.y(2:end-1); z = data.z(2:end-1);
[X,Y,Z]=ndgrid(x,y,z);
data.V = data.Potential(X,Y,Z);

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
            data.Nz = (data.zmax-data.zmin)*2^NN+1;
            data.dx = (data.xmax-data.xmin)/(data.Nx-1);
            data.dy = (data.ymax-data.ymin)/(data.Ny-1);
            data.dz = (data.zmax-data.zmin)/(data.Nz-1);
            data.x = (data.xmin:data.dx:data.xmax)';
            data.y = (data.ymin:data.dy:data.ymax)';
            data.z = (data.zmin:data.dz:data.zmax)';
            % Potential function
            x = data.x(2:end-1); y = data.y(2:end-1); z = data.z(2:end-1);
            [X,Y,Z]=ndgrid(x,y,z);
            data.V = data.Potential(X,Y,Z);
            % Initial value
            if(dN==1 && dvep==1)   
                X0=MGPE_FD3d_Init(data); 
            elseif (dN==1 && dvep>1)
                Rho_n=Rho(1:2^(NN_end-dN+1):end,1:2^(NN_end-dN+1):end,1:2^(NN_end-dN+1):end);
                X0=Rho_n(2:end-1,2:end-1,2:end-1);
            else
                X0=refine_3D(Rho);            
            end
            % update not only Rho, but also data by new values of NN
            [E,Rho] = MGPE_FD3d_Exam1(data,vep,NN,X0);
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

%% 3D plot (BOX)
% data.beta = 100; data.delta = 100; vep=1e-4; Pot='Box';
% fname_mat = strcat('MGPE_OptRegE_3D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep),'_',Pot,'.mat');
% load(fname_mat)
% x = data.x(2:end-1); y = data.y(2:end-1); z = data.z(2:end-1);
% [X,Y,Z]=ndgrid(x,y,z);
% 
% 
% fig=figure;
% sRho=smooth3(Rho,'box',5); sX=smooth3(X,'box',5);
% sY=smooth3(Y,'box',5);
% sZ=smooth3(Z,'box',5);
% h = patch(isosurface(sX,sY,sZ,sRho(2:end-1,2:end-1,2:end-1),1e-2),...
%    'FaceColor','blue','EdgeColor','none');%,...
% %    'AmbientStrength',.2,...
% %    'SpecularStrength',.7,...
% %    'DiffuseStrength',.4);
% isonormals(sRho,h)
% 
% 
% patch(isocaps(sX,sY,sZ,sRho(2:end-1,2:end-1,2:end-1),1e-2),...
%    'FaceColor','interp',...
%    'EdgeColor','none');
% colormap jet; 
% % caxis([0.01,0.03]); 
% hb=colorbar; 
% % set(hb,'YTick',[0.01:0.003:0.017]) % box1
% % set(hb,'YTick',[0.01:0.004:0.028]) % box3
% 
% title(strcat('(b) $\beta=',int2str(data.beta),',\,\delta=',int2str(data.delta),'$'),'Interpreter','latex');
% 
% daspect([1,1,1]); view(3); 
% axis tight; camlight; lighting gouraud;
% set(gca,'FontSize',25);
% xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
% zlabel('$z$','Interpreter','latex');
% 
% fname = strcat('MGPE_OptRegE_3D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep,'%1.e'),'_',Pot); % Har: harmonic potential Opt: optical potential
% saveas(fig,strcat(fname,'.fig'));save(strcat(fname,'.mat'),'E','Rho','data','vep');
% % % save .eps
% % print(fig,'-depsc',strcat(fname,'.eps')); % do not work!?
% print(fig,strcat(fname,'.png'),'-dpng');

%% 3D plot (original)
%======================================
% % comment this part if run the whole program
% data.beta = 1; data.delta = 1; vep=1e-4; Pot='Opt';
% fname_mat = strcat('MGPE_OptRegE_3D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep,'%1.e'),'_',Pot,'.mat');
% load(fname_mat)
%======================================

fig=figure;
sRho=smooth3(Rho,'box',5); iso_value=1e-2; % if there is no isovalue in fname, it is default value 1e-2
[Xf,Yf,Zf]=ndgrid(data.x,data.y,data.z);
% sX=smooth3(Xf,'box',5);
% sY=smooth3(Yf,'box',5);
% sZ=smooth3(Zf,'box',5);
% fv=isosurface(sX,sY,sZ,sRho,2e-2);
fv=isosurface(Xf,Yf,Zf,Rho,iso_value);
p1 = patch(fv,'FaceColor','green','EdgeColor','none');

daspect([1,1,1]); view(3); 
axis tight; camlight; lighting gouraud;
set(gca,'FontSize',30);
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
title(strcat('(a) $\beta=',int2str(data.beta),',\,\delta=',int2str(data.delta),'$'),'Interpreter','latex');

fname = strcat('MGPE_OptRegE_3D_',num2str(data.beta),'_',num2str(data.delta),'_',num2str(vep,'%1.e'),'_',num2str(iso_value,'%1.e'),'_',Pot); % Har: harmonic potential Opt: optical potential
saveas(fig,strcat(fname,'.fig'));save(strcat(fname,'.mat'),'E','Rho','data','vep');
% % save .eps
% print(fig,'-depsc',strcat(fname,'.eps')); % do not work!?
print(fig,strcat(fname,'.png'),'-dpng');
