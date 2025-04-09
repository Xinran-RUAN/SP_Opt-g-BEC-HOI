% test
clear; clc; 
% ============================================================
% 初始化data参数
data.beta = 1e1;
data.delta = 1e2;
data.vep = 1e-4; 
% 计算区域
Lx = 16;
data.dx = 1 / 2^3;
data.xmax = Lx;
data.xmin = -Lx;
data.Nx = (data.xmax-data.xmin) / data.dx;
data.x = (data.xmin:data.dx:data.xmax-data.dx)';
% 外势
data.Potential = @(x) 10*x.^2/2;
data.V = data.Potential(data.x);
%=================================
mass = 1;
domain_size = data.xmax - data.xmin;
x = data.x;
Xc = exp(-x.^2)/(pi^(1/2));
G = GradObj_Func(Xc, data);
L = ceil(max(abs(G)));
Xt = Xc - G / L;
Xn = myPosConProj(Xt, mass, domain_size, 1e-12);
figure(1)
plot(Xc); hold on;
plot(Xt); hold on;
plot(Xn); hold off;

figure(2)
plot(G)
