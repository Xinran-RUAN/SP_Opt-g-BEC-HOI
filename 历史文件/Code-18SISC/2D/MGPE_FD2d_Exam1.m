function [f,Rho] = MGPE_FD2d_Exam1(data,vep,NN,X0)
% opts.alg='TS';
% opts.alpha=1-1e-8; 
opts = struct('restart',-Inf,'tol',1e-8,'maxits',1e6);
% Parameters
data.vep = vep;
data.Nx = (data.xmax-data.xmin)*2^NN+1;
data.Ny = (data.ymax-data.ymin)*2^NN+1;
% Prepare data
data = MGPE_FD2d_Data(data);


% OptM solve
X0_vec=reshape(X0,(data.Nx-2)*(data.Ny-2),1);
my_MGPE_FD2d_Func = @(varargin)MGPE_FD2d_Func( data, varargin{:});
X=tfocs(my_MGPE_FD2d_Func,[],proj_simplex(2^(2*NN)),X0_vec,opts);
X_MAT=reshape(X,data.Nx-2,data.Ny-2);




% Print Energy and Chemical Potential
[f,~] = MGPE_FD2d_Func(data,X);
[M,N]=size(X_MAT);
Rho = zeros(M+2,N+2);
Rho(2:end-1,2:end-1)=X_MAT;
