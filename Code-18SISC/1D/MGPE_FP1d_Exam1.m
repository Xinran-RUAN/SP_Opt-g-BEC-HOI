function [f,Rho] = MGPE_FP1d_Exam1(data,vep,NN,X0)
% opts.alpha=1-1e-8; 
% opts = struct('restart',1000,'tol',1e-16,'alg','AT');
opts = struct('restart',-Inf,'tol',1e-20,'maxits',1e6,'alg','AT');
% Parameters
data.vep = vep;
data.Nx = (data.xmax-data.xmin)*2^NN+1;
% Prepare data
data = MGPE_FP1d_Data(data);

% OptM solve
my_MGPE_FP1d_Func = @(varargin)MGPE_FP1d_Func( data, varargin{:});
X=tfocs(my_MGPE_FP1d_Func,[],proj_my_FP((data.xmax-data.xmin)),X0,opts);

% Print Energy and Chemical Potential
[f,~] = MGPE_FP1d_Func(data,X);
Rho=[X;X(1)]; % periodic boundary condition
