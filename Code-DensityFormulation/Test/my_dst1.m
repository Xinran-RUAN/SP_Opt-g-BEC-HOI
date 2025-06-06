function f_hat = my_dst1(f) 
N = length(f) + 1;
j = (1:N-1)';     % spatial indices
l = (1:N-1);      % mode indices
S = sin(pi * j * l / N);  % (N-1 x N-1) DST-I basis

% Forward DST-I
f_hat = (2 / N) * (S' * f);  % (N-1 x 1)