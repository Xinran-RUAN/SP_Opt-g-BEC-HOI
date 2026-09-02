function [f, g] = MGPE_FP1d_Func(data,X)

beta = data.beta;
delta = data.delta;
h = data.dx;
V = data.V;
xmax=data.xmax;
xmin=data.xmin;
vep =data.vep;
N = data.Nx-1;
Lx=xmax-xmin;



%============================================================================
% kinetic energy term -- form 2
%   1/8*\frac{|\delta_x^+ \rho_j|^2}{\sqrt{\rho_{j+1/2}^2+\vep^2}}
%============================================================================
% % periodic BC (Fourier Pseudospectral)
% my_index = [ 0:N/2-1, -N/2:-1]';
my_index = [ 0:N/2-1, 0, -N/2+1:-1]'; % ?? check https://math.stackexchange.com/questions/740840/derivative-of-function-using-discrete-fourier-transform-matlab
my_lmd = 1i*2*pi/Lx*my_index;
my_xi = exp(1i*pi/N); 
my_lmd_shift = [ my_xi.^(0:N/2-1), my_xi.^(-N/2:-1)].'; % !: transpose .' not '
% 
% X_shift = real(ifft(my_lmd_shift.*fft(X))); % !!!: No ganruatee that X_shift>=0 and real
% dX_shift = real(ifft(my_lmd.*my_lmd_shift.*fft(X)));
dX = real(ifft(my_lmd.*fft(X)));
X_full = interpft(X,2*N);
dX_full = interpft(dX,2*N);
X_shift = X_full(2:2:end);
dX_shift = dX_full(2:2:end);

F1_temp = dX_shift./(sqrt(X_shift.^2+vep^2));
F2_temp = X_shift./(sqrt(X_shift.^2+vep^2));
g1 = h/8*(2*real(fft(my_lmd.*my_lmd_shift.*ifft(F1_temp)))-real(fft(my_lmd_shift.*ifft(F2_temp.*F1_temp.^2))));
d2X = real(ifft(my_lmd.^2.*fft(X)));

g = g1+h*(V.*X./sqrt(X.^2+vep^2) + X*beta - delta*d2X);
f = h/8*sum(dX_shift.*F1_temp)+h*sum(V.*(sqrt(X.^2+vep^2)-vep)+X.*X*beta/2-0.5*delta*X.*d2X);
