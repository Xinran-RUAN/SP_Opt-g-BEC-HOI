% fast computation of DCT-II(not exactly):
%   e=(e_1,e_2,...,e_N) => d=(d_1,d_2,...,d_{N-1})
% Define B=(b_jk), b_jk=cos((j-1/2)*k*pi/N)
%   j=1,2,...,N; k=1,2,...,N-1
% INPUT: 
%       e,N: [N,1]=size(e)
% OUTPUT:
%       d=B^T*e
function [d] = dct_ii(e,N)
ee=reshape([zeros(size(e)),e]',2*N,1);
E=[ee/2;0;ee(end:-1:2)/2];
dd=ifft(E)*(4*N);
d=dd(2:N);