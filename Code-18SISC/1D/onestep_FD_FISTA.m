% minimize the discrete energy by using rho
% Use APG(FISTA) method for each step
function[E,E_v,ep_test,rho]=onestep_FD_FISTA(beta,delta,N,rho,V,dh,vep,tol)
format long;
ep_test=vep;    
disp([' vep=', num2str(vep),'; dh=', num2str(dh)])

%=================================================================
% EXPREESION of ENERGY: E=0.5*rho'*H*rho+dh*V'*rho+f(rho)
%=================================================================
B=1/dh^2*(spdiags([-ones(N,1) 2*ones(N,1) -ones(N,1)], -1:1,N,N));
I=spdiags(ones(N,1), 0,N,N);
H_s=beta*dh*I+dh*delta*B;

% compute energy for given initial solution
% use Eo_v to record the energy in last step
sq_rho=sqrt(rho+vep); 
f_rho=sum((sq_rho(2:N)-sq_rho(1:N-1)).^2/dh^2)*dh/2+1/2/dh*((sq_rho(1)-sqrt(vep))^2+(sq_rho(N)-sqrt(vep))^2);
Eo_v=0.5*rho'*H_s*rho+dh*V'*rho+f_rho;


err=1;   tt=1; 
beta_step=0.8;  L_step=1; 

% initialize yy
yy=rho;
E_vep_yy=Eo_v;

while err>tol

    %==============================================================================
    % Projected Gradient Method
    %==============================================================================
    %=================================================================
    % NONLINEAR PART: f(rho)\approx f(rho^k)+grad f(rho^k)'*(rho-rho^k)
    %=================================================================



    % first choice of the kinetic energy
    sq=[sqrt(vep);sqrt(yy+vep);sqrt(vep)];  % length: N+2
    % grad (nonlinear part)
    grad_rho=-0.5/dh*diff(sq,2)./sq(2:N+1);            
    % grad E_v(yy)
    grad_E=H_s*yy+dh*V+grad_rho; 


    N_iter=1; En_vep=1;  rho_n=rho; 

    % update rho^n to rho^{n+1}
    while( N_iter==1 || Q_L<En_vep)
        if (N_iter==1) 
            L_step_n=L_step;
        end


        %=================================================================
        % CALL eplb: get t_rho(rho_n) from yy
        %=================================================================    
        grad_Er=yy-1/L_step_n*grad_E;
        rho_n = ProjectOntoL1Ball(grad_Er-min(grad_Er)+1/dh/N, 1/dh);
%         rho_n=eplb(grad_Er-min(grad_Er)+1/dh,N,1/dh,1/N);
        % compute E_v(t_tho) and Q_L(t_rho,yy)
        sq_rho_n=[sqrt(vep);sqrt(rho_n+vep);sqrt(vep)]; 
        f_rho_n=sum(diff(sq_rho_n).^2)/dh/2;
        En_vep=0.5*rho_n'*H_s*rho_n+dh*V'*rho_n+f_rho_n;
        Q_L=E_vep_yy+(rho_n-yy)'*grad_E+L_step_n/2*sum((rho_n-yy).^2);          


        % update L
        L_step_n=L_step_n/beta_step;
        N_iter=N_iter+1;
        if (L_step_n>1e16)
            fprintf('L too large! Exit loop.')
            rho_n=rho; % make error to be 0
            break
        end
    end
%     fprintf('Energy (with epsilon) is %1.10e\n',En_vep)

%     err=abs(Eo_v-En_vep)/tt; % wrong 
%     err=abs(Eo_v-En_vep); 
    err=max(abs(rho-rho_n));
    disp(err);
    
    % update E, L, yy, rho, tt
    tt_n=(1+sqrt(1+4*tt^2))/2;
    L_step=L_step_n*beta_step;
    yy=rho_n+(tt-1)/tt_n*(rho_n-rho);
    rho=rho_n;
    Eo_v=En_vep;
    tt=tt_n;



    % compute E_v(yy)
    sq_rho_yy=[sqrt(vep);sqrt(yy+vep);sqrt(vep)]; 
    f_rho_yy=sum(diff(sq_rho_yy).^2)/dh/2;
    E_vep_yy=0.5*yy'*H_s*yy+dh*V'*yy+f_rho_yy;


end

% compute E(rho) and E_v(rho)
sq_rho=sqrt(rho);
E=0.5*rho'*H_s*rho+dh*V'*rho+sum((sq_rho(2:N)-sq_rho(1:N-1)).^2/dh^2)*dh/2+1/2/dh*((sq_rho(1))^2+(sq_rho(N))^2);
fprintf('Energy is %.8f\n',E)
mu=0.5*rho'*H_s*rho+E;
fprintf('chemical potential is %.8f\n',mu)    

sq_rho=sqrt(rho+vep);
E_v=0.5*rho'*H_s*rho+dh*V'*rho+sum((sq_rho(2:N)-sq_rho(1:N-1)).^2/dh^2)*dh/2+1/2/dh*((sq_rho(1)-sqrt(vep))^2+(sq_rho(N)-sqrt(vep))^2);
fprintf('Energy(with epsilon) is %.8f\n',E_v)
mu_v=0.5*rho'*H_s*rho+E_v;
fprintf('chemical potential(with epsilon) is %.8f\n',mu_v)   



    
    


