function op = proj_my_FP( q )
op = @(varargin)proj_my_FP_detail( q, varargin{:} );

function [ v , X ] = proj_my_FP_detail( q , X , t)

%PROJ_my_SP	Projection onto the affine space 
%    X: values at the grid points
%    Aim to get the projection of q onto the space such that
%       \int X dx =q, where we treat X as a linear combination
%       of the basis sin(l\pi(x-a)/(b-a)
if nargin > 2 && t > 0
    N = length(X);
    e1 = zeros(1,N); e1(1)=1;
    mm = fft(e1)/N;
    X = X + (1/q-mm*X)/sum(mm.*mm)*mm'; % mm: row vector; X: col vector
    v=0;
else
    v = Inf;
end

% seem can be replaced by proj_simplex?? 
% b/c the coefficient = FD sum??