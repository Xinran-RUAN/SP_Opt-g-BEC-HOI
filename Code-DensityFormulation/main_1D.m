% Main file
% ============================================================ 
% 该程序以FFT离散化能量
% 使用基于梯度的一阶优化方法计算基态
% 使用投影处理可能出现的负情形
% ============================================================
clear; clc; 
tic;    % 计时开始
% ============================================================
% 初始化data参数
data.beta = 1e1;
data.delta = 1e1;
data.vep = 1e-3; 
% 计算区域
Lx = 32;
data.dx = 1/2^5;
data.xmax = Lx;
data.xmin = -Lx;
data.Nx = (data.xmax-data.xmin) / data.dx;
data.x = (data.xmin:data.dx:data.xmax-data.dx)';
% 外势
data.Potential = @(x) x.^2/2;
data.V = data.Potential(data.x);

% ============================================================
% 初始化数值解
x = data.x;
X0 = exp(-x.^2)/(pi^(1/2));

% ============================================================
% 计算基态解与基态能量
[E,Rho] = Solve_Ground_State(X0, data);           

% ============================================================
toc;    % 计时结束
% ============================================================
% 画图
plot(data.x,Rho)

%% save data
folderPath = './Data';  % 可以是绝对路径或相对路径

if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end

cd(folderPath);

% filename = strcat('MGPE-FP1d-Bet-',int2str(data.beta),...
%     '-Del-',int2str(data.delta),...
%     '-Vep-',num2str(data.vep),...
%     '-dx-',num2str(data.dx),'.mat');
filename = strcat('MGPE-FP1d-Bet-',int2str(data.beta),...
    '-Del-',int2str(data.delta),...
    '-Vep-',num2str(data.vep),...
    '-dx-',num2str(data.dx),'_reg1_FISTA_CD_highE.mat'); % a new regularization of the kinetic part
save(filename)
cd ..
