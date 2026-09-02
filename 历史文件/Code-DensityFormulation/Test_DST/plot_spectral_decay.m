function plot_spectral_decay(rho)
    % rho: 输入的是一个实值向量，代表离散函数值，长度 N
    N = length(rho);
    
    % 计算傅里叶变换（不做缩放）
    hat_rho = fft(rho);
    
    % 将频率居中排列：低频在中间
    hat_rho_shifted = fftshift(hat_rho);

    % 计算模长（幅度谱）
    abs_hat = abs(hat_rho_shifted);

    % 构造频率轴
    k = -N/2:N/2-1;  % 假设N为偶数

    % 绘图（对数坐标）
    figure;
    semilogy(k, abs_hat, 'b-o', 'LineWidth', 1.5);
    xlabel('Wave number k');
    ylabel('|\hat{\rho}_k| (log scale)');
    title('Spectral Coefficient Decay');
    grid on;
end