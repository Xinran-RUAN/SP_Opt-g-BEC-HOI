% PLOT
%% Fix beta 
% beta=10; delta=[10;40;500];
% vep=0.1./2.^(0:5)';
% 
% err1=[0.025903044058078
%    0.020745508423675
%    0.015648536115150
%    0.011161544002294
%    0.007575885335606
%    0.004929196075990];
% 
% err2=[0.016193165928837
%    0.013362876073820
%    0.010440550906344
%    0.007740973726133
%    0.005465701671859
%    0.003690432631236];
% 
% err3=[0.003919299924070
%    0.003390375193355
%    0.002796683047233
%    0.002205429062721
%    0.001669810087535
%    0.001218354706200];
% 
% fig = figure(1);
% Linewidth=3; Markersize=17;
% set(fig, 'Position', [360 278 840 560])
% loglog(vep,err1,'-d','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold on;
% loglog(vep,err2,'-.s','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold on;
% loglog(vep,err3,'--o','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold off;
% 
% textsize=30;
% % axis tight;
% set(gca,'FontSize',textsize,'Linewidth',2);
% set(gca,'XTick',[vep(1)/32,vep(1)/16,vep(1)/8,vep(1)/4,vep(1)/2,vep(1)]);
% set(gca,'XTickLabel',{'$\varepsilon$/32';'$\varepsilon$/16';'$\varepsilon$/8';'$\varepsilon$/4';'$\varepsilon/2$';'$\varepsilon=0.1$'},'TickLabelInterpreter','latex');
% set(gca,'YTick',[0.001,0.002,0.005,0.01,0.02]);
% set(gca,'YTickLabel',{'$0.1\%$';'$0.2\%$';'$0.5\%$';'$1\%$';'$2\%$'},'TickLabelInterpreter','latex');
% xlim([0.0027,0.11]); 
% ylim([0.001, 0.03]);
%  title('$\beta=10$','Interpreter','latex')
% 
% % xlabel('$\varepsilon$','Interpreter','latex','FontSize',textsize)
% ylabel('$\frac{|E^{\varepsilon}-E_{ex}|}{E_{ex}}$','Interpreter','latex','FontSize',textsize)
% 
% legend_label={'$\delta=10$','$\delta=40$','$\delta=100$'};
% legend(legend_label,'Location','southeast','Interpreter','latex')
% 
% fname = 'relative_error1';
% saveas(fig,strcat(fname,'.fig'));
% % % save .eps
% % print(fig,'-depsc',strcat(fname,'.eps'));

%% Fix delta 
vep=[1e-2,1e-3,1e-4];
delta=10; beta=[10;20;40;80;160;320];

err1=[0.025903044058078
   0.020745508423675
   0.015648536115150
   0.011161544002294
   0.007575885335606
   0.004929196075990];

% err2=[0.008580993697359
%    0.007156166714299
%    0.005620745560681
%    0.004162688524254
%    0.002921198017842
%    0.001956047061340];
% 
% err3=[0.003436221577158
%    0.002969235622122
%    0.002423489422286
%    0.001864057313836
%    0.001354423404034
%    0.000934994696826];

fig = figure(1);
Linewidth=3; Markersize=17;
set(fig, 'Position', [360 278 840 560])
loglog(beta,err1,'-d','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold on;
% loglog(vep,err2,'-.s','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold on;
% loglog(vep,err3,'--o','Linewidth',Linewidth,'MarkerSize',Markersize);
% hold off;

textsize=30;
% axis tight;
set(gca,'FontSize',textsize,'Linewidth',2);
set(gca,'XTick',[beta(1),beta(2),beta(3),beta(4),beta(5),beta(6)]);
% set(gca,'XTickLabel',{'$\varepsilon$/32';'$\varepsilon$/16';'$\varepsilon$/8';'$\varepsilon$/4';'$\varepsilon/2$';'$\varepsilon=0.1$'},'TickLabelInterpreter','latex');
% set(gca,'YTick',[0.001,0.002,0.005,0.01,0.02]);
% set(gca,'YTickLabel',{'$0.1\%$';'$0.2\%$';'$0.5\%$';'$1\%$';'$2\%$'},'TickLabelInterpreter','latex');
% xlim([0.0027,0.11]); 
% ylim([0.0008, 0.03]);
 title('$\delta=10$','Interpreter','latex')

% xlabel('$\varepsilon$','Interpreter','latex','FontSize',textsize)
ylabel('$\frac{|E^{\varepsilon}-E_{ex}|}{E_{ex}}$','Interpreter','latex','FontSize',textsize)

legend_label={'$\varepsilon=10^{-4}$'};
legend(legend_label,'Location','southeast','Interpreter','latex')

fname = 'relative_error2';
saveas(fig,strcat(fname,'.fig'));
% % save .eps
% print(fig,'-depsc',strcat(fname,'.eps'));