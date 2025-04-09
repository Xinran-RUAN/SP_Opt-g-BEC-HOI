function [E, X] = my_tfocsAT(X0, data)
% OptM solve
my_MGPE_FP1d_Func = @(varargin)myDiscreteEnergy(data, varargin{:});
opts = struct('restart',-Inf,'tol',1e-16,'alg','AT');
X=tfocs(my_MGPE_FP1d_Func,[],myPosConProj(?), X0, opts);

% Print Energy and Chemical Potential
[E,~] = my_discrete_energy(data,X);

