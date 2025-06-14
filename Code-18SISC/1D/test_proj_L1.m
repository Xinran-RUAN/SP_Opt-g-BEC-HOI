clc; clear;
N=99; x=((1:N)-50).^2'; fun=@fun_test;
%-------------------------------------------------------------------------
% Optimization:
% Model :
%   min f(X), s.t., ||X||_1 = radius, X>=0, where X \in R^{n}
%       g(X) = grad f(X)
% Output :
%   x --- solution
% Test example: f(X)=sum(X_j^2) s.t. sum(X_j)=1, X_j>=0
% g(X)=2*X
% parameters for control
opts.tol = 1e-15;      
opts.radius = 1;      
opts.maxit = 50000;      
%-------------------------------------------------------------------------
% copy parameters
tol = opts.tol;
radius = opts.radius;
maxit = opts.maxit;
%% main iteration
tt=1; N=length(x);
[E_yy,grad_E]=feval(fun, x);
h=1;
% initialize yy & rho
rho=x; yy=x; 
L_step=100; beta_step=0.999; % back-tracking
for iter = 1 : maxit 
    %==============================================================================
    % Projected Gradient Method
    %==============================================================================
     En=1;   
    %===========================================
    % update rho^n to rho^{n+1} (one update)
    % determine L using backtracking
    %===========================================
    N_iter=1;rho_n=rho;
    while( N_iter==1 || Q_L<En)
        if (N_iter==1) 
            L_step_n=L_step;
        end
        %=================================================================
        % CALL Project: get t_rho(rho_n) from yy
        %=================================================================    
        grad_Er=yy-1/L_step_n*grad_E;
        rho_n = ProjectOntoL1Ball(grad_Er-min(grad_Er)+radius/h/N, radius/h);
        %===========================================
        % compute E_v(t_tho) and Q_L(t_rho,yy)
        %===========================================
        [En,~] = feval(fun, rho_n);
        Q_L=E_yy+(rho_n-yy)'*grad_E+L_step_n/2*sum((rho_n-yy).^2); 
        %===========================================
        % update L
        %===========================================
        L_step_n=L_step_n/beta_step;
        N_iter=N_iter+1;      
    end
    
    L_step=L_step_n*beta_step;
    
    %==============================
    %   termination
    %==============================
    err=max(abs(rho-rho_n));
    if  err< tol
        break;
    end
    disp(err)
    %==============================
    % update E, L, yy, rho, tt
    %==============================
    tt_n=(1+sqrt(1+4*tt^2))/2;
    yy=rho_n+(tt-1)/tt_n*(rho-rho_n);
    rho=rho_n;
    tt=tt_n;
    [E_yy,grad_E] = feval(fun, yy);
end
   

plot(1:N,rho)

 