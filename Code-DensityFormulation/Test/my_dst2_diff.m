%% dst2: 
% 输出：半点处的导数值
function [dX_h] = my_dst2_diff(X, L)
X = X(:); % ensure column vector
N = length(X) + 1;  % Total N grid intervals, N-1 interior points

% Create sin basis matrix for DST-I
j = (1:N-1)';     % spatial indices
l = (1:N-1);      % mode indices
S = sin(pi * j * l / N);  % (N-1 x N-1) DST-I basis

% Forward DST-I
X_hat = (2 / N) * (S' * X);  % (N-1 x 1)

% Multiply by wave number for derivative
nu = (pi / L) * (1:N-1);  % (1 x N-1)
dX_hat = nu(:) .* X_hat;

% Derivative basis is cos(pi * l * j_h / N)
j_h = (1/2:N-1/2)';
C = cos(pi * j_h * l / N);  % (N x N-1)
dX_h = C * dX_hat;      % (N-1 x 1)
end
