% fast computation of DST-III(not exactly):
%   e=(e_1,e_2,...,e_N) => d=(d_1,d_2,...,d_{N-1}) 
% Define D=(d_jk), d_jk=sin((j-1/2)*k*pi/N)
%   j=1,2,...,N; k=1,2,...,N-1
% INPUT: 
%       d,N: [N,1]=size(e)
% OUTPUT:
%       e=D*d
function [e] = dst_iii(d,N)
dd=[0;d];
D=[dd/2i;zeros(2*N+1,1);-dd(end:-1:2)/2i];
ee=ifft(D)*(4*N);
e=ee(2:2:2*N);