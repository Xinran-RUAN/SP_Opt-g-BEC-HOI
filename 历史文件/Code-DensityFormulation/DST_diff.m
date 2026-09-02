function [dX] = DST_diff(X, L)
%%
    N = length(X);
    mu = pi / L * (1:N);
    mu = reshape(mu, size(X)); % Lambda与X形状保持一致
%%
    X_hat = dst(X);
    dX_hat = mu .* X_hat;
    dX = idct(dX_hat);
end