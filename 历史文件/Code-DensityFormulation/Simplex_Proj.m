function[y_out] = Simplex_Proj(x_in, mass)
x = x_in / mass; 
    % 将向量排序
    [u, ~] = sort(x, 'descend');
    % 计算累积和
    cumsum_u = cumsum(u);
    % 找到满足条件的最大的j
    rho = find(u > (cumsum_u - 1) ./ (1:length(u))', 1, 'last');
    % 计算lambda值
    lambda = (cumsum_u(rho) - 1) / rho;
    % 计算投影后的向量
    y = max(x - lambda, 0);
y_out = y * mass;
end
