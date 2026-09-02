function [f,Rho] = MGPE_FD3d_Exam1(data,vep,NN,X0)

% opts.alg='TS';
% opts.alpha=1-1e-8; 
opts = struct('restart',-Inf,'tol',1e-8,'maxits',1e6);
% Parameters
data.vep = vep;
data.Nx = (data.xmax-data.xmin)*2^NN+1;
data.Ny = (data.ymax-data.ymin)*2^NN+1;
data.Nz = (data.zmax-data.zmin)*2^NN+1;
dim_x=(data.Nx-2)*(data.Ny-2)*(data.Nz-2);
X_INIT=reshape(X0,dim_x,1);

% Prepare data
data = MGPE_FD3d_Data(data);

% OptM solve
my_MGPE_FD3d_Func = @(varargin)MGPE_FD3d_Func( data, varargin{:});
X=tfocs(my_MGPE_FD3d_Func,[],proj_simplex(2^(3*NN)),X_INIT,opts);
X_MAT=reshape(X,data.Nx-2,data.Ny-2,data.Nz-2);


% Print Energy and Chemical Potential
[f,~] = MGPE_FD3d_Func(data,X);
Rho = zeros(data.Nx,data.Ny,data.Nz);
Rho(2:end-1,2:end-1,2:end-1)=X_MAT;
