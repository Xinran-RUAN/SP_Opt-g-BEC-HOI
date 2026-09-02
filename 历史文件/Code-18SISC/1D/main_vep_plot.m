clc;clear;
beta=10; delta=10; NN=6; Linewidth=2;textsize=25;
mvep_list={1e-2/4,1e-2/2,1e-2};
mvep_name={'\frac{1}{400}','\frac{1}{200}','\frac{1}{100}'};
% Linestyle={'-bs','-rd','-mo'};
Linestyle={'--','-.','.'};
% load exact solution
filename=strcat('GPE-SP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-NN-6');
load(filename)
Rho_ex=Phi.^2; x_ex=data.x;
figure(1); 
% plot(data.x,Rho_ex,'-r','Linewidth',Linewidth); hold on;
% legend_label={'$\rho_{g}\,(\varepsilon=0)$'};
legend_label={};

for jj=1:length(mvep_list)
% for jj=length(mvep_list):-1:1
    vep=mvep_list{jj};
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep),'-NN-',int2str(NN),'_n.mat');
%     disp(filename)
    load(filename)
    rho_ex=interp1(x_ex,Rho_ex,data.x);
    plot(data.x,Rho-rho_ex,Linestyle{jj},'Linewidth',Linewidth); hold on;
%     plot(data.x,Rho-Rho_ex,Linestyle{jj},'Linewidth',Linewidth); hold on;
    legend_label=[legend_label,strcat('$\varepsilon=',mvep_name{jj},'$')];
end
% plot(data.x,Rho_ex,Linestyle{jj},'Linewidth',Linewidth);
hold off
xlim([-5,5]); legend(legend_label,'Interpreter','latex','Location','best');
% xlim([3,3.8]); ylim([0 0.035]); % subplot



set(gca,'FontSize',textsize,'Linewidth',2);
xlabel('$x$','Interpreter','latex'); 
ylabel('$\rho_{g}^\varepsilon-\rho_{g}$','Interpreter','latex');
% ylabel('$\rho_{g}^\varepsilon$','Interpreter','latex');
