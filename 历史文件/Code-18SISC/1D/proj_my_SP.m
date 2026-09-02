function op = proj_my_SP( q )
op = @(varargin)proj_my_SP_detail( q, varargin{:} );

function [ v , X ] = proj_my_SP_detail( q , X , t)

%PROJ_my_SP	Projection onto the affine space 
%    X: values at the grid points
%    Aim to get the projection of q onto the space such that
%       \int X dx =q, where we treat X as a linear combination
%       of the basis sin(l\pi(x-a)/(b-a)
if nargin > 2 && t > 0
    N = length(X)+1;
    ll = (1-(-1).^(1:N-1))./(1:N-1);
    mm = idst(ll);
    X = X + (1/q-mm*X)/sum(mm.*mm)*mm'; % mm: row vector; X: col vector
    v=0;
else
    v = Inf;
end
