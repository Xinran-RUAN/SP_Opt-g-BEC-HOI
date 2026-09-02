function X0 = MGPE_FD3d_Init(data)
init=data.init; 
if strcmp(init,'Auto')
    if (data.beta < 30 && data.delta<30)
        init = 'Gaussian';
    elseif (data.delta<30)
        init = 'TFA'; % large beta, small delta
    else
        init='TFA_D'; % large delta, ? small beta
    end
end

x = data.x(2:end-1); y = data.y(2:end-1); z = data.z(2:end-1);
hx=data.dx; hy=data.dy; hz=data.dz;
[X,Y,Z]=ndgrid(x,y,z);

% Gaussian Approximation
if strcmp(init,'Gaussian') || (strcmp(init,'TFA') && data.beta==0)
    X0 = exp(-X.^2-Y.^2-Z.^2);
% % Thomas-Fermi Approximation % only for \gm=1!
% elseif strcmp(init,'TFA') && data.beta~=0
%     beta = data.beta;
%     V = data.Potential(X,Y);
%     X0 = max(((beta*4/pi)^(1/2) - X.^2-Y.^2)/beta/2,0);
% % Thomas-Fermi Approximation(large delta)
% elseif strcmp(init,'TFA_D') && data.delta~=0
%     delta = data.delta;
%     gm=1;
%     R=(4^2*6*delta/pi/gm^2)^(1/6);
%     X0=max(gm^2*max(R^2-X.^2-Y.^2,0).^2/32/delta,0);
% 
% % Random Approximation
% else
%     X0 = rand(data.Nx-2,1);
end
X0 = X0/sum(sum(sum(X0)))/hx/hy/hz;