function [f,Rho] = MGPE_SP1d_Exam1(data,vep,NN,X0)
% opts.alpha=1-1e-8; 
opts = struct('restart',-Inf,'tol',1e-16,'alg','AT');
% opts = struct('restart',-Inf,'tol',1e-16,'maxits',1e6,'alg','AT');
% Parameters
data.vep = vep;
data.Nx = (data.xmax-data.xmin)*2^NN+1;
% Prepare data
data = MGPE_SP1d_Data(data);

% OptM solve
my_MGPE_SP1d_Func = @(varargin)MGPE_SP1d_Func( data, varargin{:});
% X=tfocs(my_MGPE_SP1d_Func,[],proj_simplex(2^NN),X0,opts); % projection?? 
%     % what is the feasible set for the SP discretization??
X=tfocs(my_MGPE_SP1d_Func,[],proj_my_SP((data.xmax-data.xmin)/pi),X0,opts);

% Print Energy and Chemical Potential
[f,~] = MGPE_SP1d_Func(data,X);
Rho=[0;X;0];
