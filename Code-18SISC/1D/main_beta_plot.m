clc;clear;
beta_list=[10,100,1000]; 
delta=10; NN_plot=6; Linewidth=2;textsize=25;
vep=1e-4;
mbeta_name={'10','100','1000'};
% Linestyle={'-bs','-rd','-mo'};
Linestyle={'--','-.','.'};

figure(1); 
% plot(data.x,Rho_ex,'-r','Linewidth',Linewidth); hold on;
% legend_label={'$\rho_{g}\,(\varepsilon=0)$'};
legend_label={};

for jj=1:length(beta_list)
    beta=beta_list(jj);
    % load exact solution
    filename=strcat('GPE-SP1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-NN-6_n');
    load(filename)
    Rho_ex=Phi.^2; x_ex=data.x;
    filename=strcat('MGPE-FD1d-Bet-',int2str(beta),'-Del-',int2str(delta),'-Vep-',num2str(vep),'-NN-',int2str(NN_plot),'_n.mat');
    load(filename)
    rho_ex=interp1(x_ex,Rho_ex,data.x,'linear','extrap');
    plot(data.x,Rho-rho_ex,Linestyle{jj},'Linewidth',Linewidth); hold on;
%     plot(data.x,Rho,Linestyle{jj},'Linewidth',Linewidth); hold on;
    % plot(data.x,Rho_ex,Linestyle{jj},'Linewidth',Linewidth); hold on;
    legend_label=[legend_label,strcat('$\beta=',mbeta_name{jj},'$')];
end
hold off
xlim([-50,50]); legend(legend_label,'Interpreter','latex','Location','best'); 
% ylim([-1e-2 0.3]); 



set(gca,'FontSize',textsize,'Linewidth',2);
xlabel('$x$','Interpreter','latex'); 
ylabel('$\rho_{g}^\varepsilon-\rho_{g}$','Interpreter','latex'); set(gca,'xtick',[-50:20:50]);
% ylabel('$\rho_{g}^\varepsilon$','Interpreter','latex');

annotation(figure(1),'textbox',...
    [0.184928571428571 0.771428571428574 0.103529575892857 0.0785714285714286],... 
    'String','(c)',...%     'String','(a)',...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica');

