function [rho, grad_E]= proj_L1(x, fun,opts, data)
%-------------------------------------------------------------------------
% Optimization:
%
% Model :
%   min f(X), s.t., ||X||_1 = radius, X>=0, where X \in R^{n}
%       g(X) = grad f(X)
% Input:
%           X --- ||X||_2 = radius
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                 
%        opts --- option structure with fields:
%                 record = 0, no print out
% Output:
%           x --- solution
%           g --- gradient of x
%
% Change epph.h for the following stop criteria:
%                 mxitr       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
% -------------------------------------
% -------------------------------------
% Author: Xinran Ruan
%   adapted from OptManiMulitBallGBB.m by Zaiwen Wen
%-------------------------------------------------------------------------
% parameters for control
if ~isfield(opts, 'tol');    opts.tol = 1e-10;      end
if ~isfield(opts, 'radius');    opts.radius = 1;      end
if ~isfield(opts, 'maxit');    opts.maxit = 50000;      end
%-------------------------------------------------------------------------------
% copy parameters
tol = opts.tol;
radius = opts.radius;
maxit = opts.maxit;
%% main iteration
tt=1; N=length(x);
[E_yy,grad_E]=feval(fun, x, data);
h=data.dx; Eo=E_yy;
% initialize yy & rho
rho=x;
yy=x; 

L_step=100; beta_step=0.99; % back-tracking
% vep=data.vep; beta=data.beta; delta=data.delta; % no back-tracking
% L_step=4/h/vep*sqrt(1+1/h/vep)+h*beta+8*delta/h; % no back-tracking
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
    while 1
        if (N_iter==1) 
            L_step_n=L_step;
        end
        %=================================================================
        % CALL Projection: get t_rho(rho_n) from yy
        %=================================================================    
        grad_Er=yy-1/L_step_n*grad_E;
%         rho_n=eplb(grad_Er-min(grad_Er)+radius/h/N,N,radius/h,1/N);
        rho_n = ProjectOntoL1Ball(grad_Er-min(grad_Er)+2*radius/h/N, radius/h);
        if (max(abs(yy-rho_n))<1e-16)
            fprintf('yy=rho_n\n')
%             L_step_n=L_step;
            break;
        end
        %===========================================
        % compute E_v(t_tho) and Q_L(t_rho,yy)
        %===========================================
        [En,~,~] = feval(fun, rho_n, data);
        Q_L=E_yy+(rho_n-yy)'*grad_E+L_step_n/2*sum((rho_n-yy).^2); 
        if Q_L>En
            break;
        end
        %===========================================
        % update L
        %===========================================
        L_step_n=L_step_n/beta_step;
        N_iter=N_iter+1;
%         if (L_step_n>1e16)
%             fprintf('L too large! Exit loop.')
%             rho_n=rho; % make error to be 0
%             break
%         end        
    end
%     disp(iter)
    

    L_step=L_step_n*beta_step;
%     %===========================================
%     % update rho^n to rho^{n+1} (one update)
%     % determine L without backtracking
%     %===========================================  
%     grad_Er=yy-1/L_step*grad_E;
%     rho_n=eplb(grad_Er-min(grad_Er)+radius/h,N,radius/h,1/N);
    
    %==============================
    %   termination
    %==============================
    err=sqrt(sum(h*(rho-rho_n).^2));
    if err<1.2E-8
        pause(0.1)
    end
    if iter==maxit
        fprintf('iteration = maxit\n')
    end
%     err=abs(Q_L-En);
    if  err< tol
%         Eo-En
        break;
    end
%     disp(err)
    %==============================
    % update E, L, yy, rho, tt
    %==============================
    tt_n=(1+sqrt(1+4*tt^2))/2;
    yy=max(rho_n+(tt-1)/tt_n*(rho_n-rho),0);
    rho=rho_n; Eo=En; 
    tt=tt_n;
    [E_yy,grad_E,~] = feval(fun, yy, data);
end
   

X=rho; h=data.dx; V=data.V; beta=data.beta; delta=data.delta; D=data.D;
sq=[0;sqrt(X);0];
E_wo_vep = h*sum(V.*X+X.*X*beta/2+delta/2*X.*(D*X))+sum(diff(sq).^2)/h/2;

% disp(iter); % disp(L_step);

fprintf('Energy(with epsilon) is %.8f, Energy(without epsilon) is %.8f\n',En, E_wo_vep)





 