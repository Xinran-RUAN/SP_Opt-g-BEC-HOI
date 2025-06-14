% fast computation of DCT-III(not exactly):
%   d=(d_1,d_2,...,d_{N-1}) => e=(e_1,e_2,...,e_N) 
% Define B=(b_jk), b_jk=cos((j-1/2)*k*pi/N)
%   j=1,2,...,N; k=1,2,...,N-1
% INPUT: 
%       d,N: [N,1]=size(e)
% OUTPUT:
%       e=B*d
function [e] = dct_iii(d,N)
dd=[0;d];
D=[dd/2;zeros(2*N+1,1);dd(end:-1:2)/2];
ee=ifft(D)*(4*N);
e=ee(2:2:2*N);