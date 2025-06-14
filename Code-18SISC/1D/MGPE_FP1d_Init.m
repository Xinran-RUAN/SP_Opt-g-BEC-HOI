function X0 = MGPE_FP1d_Init(data)
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

x = data.x(1:end-1);
h=data.dx;

% Gaussian Approximation
if strcmp(init,'Gaussian') || (strcmp(init,'TFA') && data.beta==0)
    X0 = exp(-x.^2)/(pi^(1/2));
% Thomas-Fermi Approximation
elseif strcmp(init,'TFA') && data.beta~=0
    beta = data.beta;
    V = data.Potential(x);
    X0 = max((1/2*(beta*3/2)^(2/3) - V)/beta,0);
% Thomas-Fermi Approximation(large delta)
elseif strcmp(init,'TFA_D') && data.delta~=0
    delta = data.delta;
    gm=1;
    R=(45*delta/2/gm^2)^(1/5);
    X0=max(gm^2*max(R^2-x.^2,0).^2/24/delta,0);
% % Excited 1
% elseif strcmp(init,'Excited1')
%     % Phi = 1/pi^(d/4)*gamma_x^(1/4)
%     %       * exp(-(gamma_x*x^2)/2)
%     %       * sqrt(2)*x
%     X0 = exp(-x.^2/2)/(pi^(1/4)) ...
%          .* (sqrt(2)*x);

% Random Approximation
else
    X0 = rand(data.Nx-2,1);
end
X0 = X0/norm(X0,1)*1/h;
