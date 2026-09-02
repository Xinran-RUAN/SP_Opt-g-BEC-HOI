function [dX_T] = DST_diff_T(X, L)
    %% cut -- transpose
    [m, n] = size(X);
    if n == 1 
        X_ext = [zeros(m+2, 1); X]; % N+3至2N+2
    else
        X_ext = [zeros(1, n+2), X];
    end
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
    dX_ext_T = fourier_diff_T(X_ext, Lambda);
    %%  odd extension -- transpose
    dX_T = OddExtension_T(dX_ext_T);
end