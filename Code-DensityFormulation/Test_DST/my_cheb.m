function [D, x] = my_cheb(N)
% Compute Chebyshev differentiation matrix D and grid points x
% N: number of intervals → N+1 grid points

if N==0
    D=0; x=1;
    return
end

x = cos(pi * (0:N) / N)';  % Chebyshev points, column vector
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';  % weights

X = repmat(x,1,N+1);
dX = X - X';

D = (c * (1./c)')./(dX + (eye(N+1))); % Off-diagonal entries
D = D - diag(sum(D,2));               % Diagonal entries
end
