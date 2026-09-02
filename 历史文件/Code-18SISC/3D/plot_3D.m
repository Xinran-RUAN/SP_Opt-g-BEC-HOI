% 去除边界点
x = data.x(2:end-1); y = data.y(2:end-1); z = data.z(2:end-1);
[X,Y,Z]=ndgrid(x,y,z);

% 平滑处理数据
fig=figure;
sRho=smooth3(Rho,'box',5); sX=smooth3(X,'box',5);
sY=smooth3(Y,'box',5);
sZ=smooth3(Z,'box',5);

% 绘制等值面
h = patch(isosurface(sX,sY,sZ,sRho(2:end-1,2:end-1,2:end-1),1e-2),...
   'FaceColor','blue','EdgeColor','none');
isonormals(sRho,h)

% 绘制截面
patch(isocaps(sX,sY,sZ,sRho(2:end-1,2:end-1,2:end-1),1e-2),...
   'FaceColor','interp',...
   'EdgeColor','none');
colormap jet; 
hb=colorbar; 

% 标题、坐标等设置
title(strcat('(b) $\beta=',int2str(data.beta),',\,\delta=',int2str(data.delta),'$'),'Interpreter','latex');

daspect([1,1,1]); view(3); 
axis tight; camlight; lighting gouraud;
set(gca,'FontSize',25);
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');