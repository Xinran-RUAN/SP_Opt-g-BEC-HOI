% n为偶数，即向量长为偶数
function [g]=myifft(Ff)
[m,n]=size(Ff);
% inverse Fourier transform
if (m-1)*(n-1)~=0
    disp('Error')
end
if m==1
    f=[Ff(n/2+1:n),Ff(1:n/2)];
    g=ifft(f)*n;
else
    f=[Ff(m/2+1:m);Ff(1:m/2)];
    g=ifft(f)*m;
end
end