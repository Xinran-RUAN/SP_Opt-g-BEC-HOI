%% 
% 奇延拓（odd extension）虽然能实现快速正弦变换，但确实可能导致精度下降，尤其是在处理导数时。
function [dX] = DST_diff(X, L)
%%  odd extension
    X_OddExt = OddExtension(X);

%%
    [m,n] = size(X);
    if n == 1
        N = m;
        Lambda = 2 * pi * 1i / (2 * L) * [(0:N)'; (-N-1:-1)'];
    else
        N = n;
        Lambda = 2 * pi * 1i / (2 * L) * [0:N, -N-1:-1];
    end

%%
    dX_OddExt = fourier_diff(X_OddExt, Lambda);
    dX = dX_OddExt(N+1:end); % N+1至2N+2: 比X多两个边界值，边界速度非0
end