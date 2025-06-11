% test DST
dx = 1;
x = (0:dx/2^2:4)';
x_in = x(2:end-1);
x_h = 0.5 * (x(1:end-1) + x(2:end));
%% test my_dst
% y = sin(4 * pi * x);
% y = x.^2 .* (4-x).^2;
y = x .* (4-x);
y_hat = dst(y(2:end-1)) * (2/(length(x)-1));
my_y_hat0 = my_dst(y(2:end-1));
my_y_hat1 = dst1_scalar(y(2:end-1));
my_y_hat2 = my_dst1(y(2:end-1));
disp('Test DST_Func')
disp(max(abs(my_y_hat1 - y_hat)));
disp(max(abs(my_y_hat1 - my_y_hat0)));
disp(max(abs(my_y_hat1 - my_y_hat2)));
% 结论：my_dst有明显误差，这表明奇延拓影响精度;对于y = 4 - x.^2，计算有误？？？

%% test y_hat
disp('Test y_hat')
N = length(x) - 2;
MAT_Sin = sin(pi / L * (x - min(x)) * (1:N));
y0 = MAT_Sin * my_y_hat0;
y1 = MAT_Sin * my_y_hat1;
y2 = MAT_Sin * my_y_hat2;
disp(max(abs(y0 - y)));
disp(max(abs(y1 - y)));
disp(max(abs(y2 - y)));
%% test DST_diff
% % 测试例子1
% % 该测试例子的一阶导数不符合Neumann边界，无法用cos列表示，因此无法得到谱精度
y = x .* (4-x);
dy_ex = 4 - 2 * x;
dy_in_ex = 4 - 2 * x_in;
dy_h_ex = 4 - 2 * x_h; 
% % 测试例子2
% % 即使是Neumann边界，高阶导数非Neumann边界也会导致无谱精度？？
% y = x.^2 .* (4-x).^2;
% dy_ex = 2 * x .* (4-x).^2 - 2 * x.^2 .* (4-x);
% dy_in_ex = 2 * x_in .* (4-x_in).^2 - 2 * x_in.^2 .* (4-x_in);
% dy_h_ex = 2 * x_h .* (4-x_h).^2 - 2 * x_h.^2 .* (4-x_h); 
% % 测试例子3
% y = sin(4 * pi * x);
% dy_ex = 4 * pi * cos(4 * pi * x);
% dy_in_ex = 4 * pi * cos(4 * pi * x(2:end-1));
% dy_h_ex = 4 * pi * cos(4 * pi * x_h);

%% 数值测试
N = length(y(1:end-1));
L = max(x) - min(x);
dy = DST_diff(y(2:end-1), L);
dy1 = my_dst1_diff(y(2:end-1), L);
dy2 = my_dst2_diff(y(2:end-1), L); % half grid

Lambda = 2 * pi * 1i / L * [(0:N/2)'; (-N/2+1:-1)'];
dy_fft = fourier_diff(y(1:end-1), Lambda);

disp('Test diff_y')
disp(max(abs(dy - dy_in_ex)));
disp(max(abs(dy1 - dy_in_ex)));
disp(max(abs(dy2 - dy_h_ex)));
disp(max(abs(dy_fft - dy_ex(1:end-1))));

plot(dy_fft - dy_ex(1:end-1)) % 误差出现在边界处

% figure(1)
% plot(x, dy, 'b-'); hold on;
% plot(x(2:end-1), dy_ex, 'r--'); hold off;


