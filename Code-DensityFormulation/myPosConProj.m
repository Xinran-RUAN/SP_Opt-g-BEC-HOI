% Compute the positive and conservative projection by semi-smooth method
function [y_out] = myPosConProj(x_in,mass_in,domain_in,myeps)
% y_out     : output function at collocation points
% x_in      : input function at collocation points, typically P_N f (x_j)
% mass_in   : desired conserved mass, which is equal to 1 by default
% domain_in : size of domain, for normalization
% myeps     : absolute convergent criterion

if nargin==1
    mass_in     = 1;
    domain_in   = 1; 
    myeps       = 1e-14;
elseif nargin==2
    domain_in   = 1; 
    myeps       = 1e-14;
elseif nargin==3
    myeps       = 1e-14;
end

NofDof = length(x_in);

x_in  = reshape(x_in, NofDof, 1);

% matrix for conservation
C_mat   = domain_in / NofDof * ones(1,NofDof); % here we require C_mat is real
C_mat_T = C_mat.'; 

% initial value for lambda
mylambda = 0;

% main loop : Newton iteration; set max iteration number be 50
for it = 1 : 50
    % Newton's method: Jacobi * delta = - Residual ===> delta = ...  

    f_plus_cTlambda = x_in + C_mat_T * mylambda;

    project_pos_vec = (f_plus_cTlambda > 0);

    Newton_Residual = mass_in - C_mat(:,project_pos_vec) * f_plus_cTlambda(project_pos_vec);
    Newton_Residual_norm = norm(Newton_Residual);

    % disp([num2str(it-1),' Newton iteration, residual = ',num2str(Newton_Residual_norm)]);
    if Newton_Residual_norm<myeps
        break;
    end

    % the below code should add regularization if it is = 0
    Newton_Jacobi = C_mat(:,project_pos_vec) * C_mat_T(project_pos_vec,:);

    if ( abs(Newton_Jacobi)/Newton_Residual_norm < myeps) 
        Newton_Jacobi = myeps * Newton_Residual_norm;
    end

    % line search if required; currently no line search
    mylambda = mylambda + Newton_Jacobi \ Newton_Residual;

end

% reconstruction of solution
y_out = x_in + C_mat_T * mylambda;
y_out(y_out<0) = 0;

end