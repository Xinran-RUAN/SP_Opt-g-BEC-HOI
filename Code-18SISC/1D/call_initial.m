% This subprogram builds the initial value of rho
function[rho]=call_initial(x,pot,initialcase,beta,delta,gm)
if (strcmp(pot,'har')&& strcmp(initialcase,'ld'))
    R=(45*delta/2/gm^2)^(1/5);
    rho=gm^2*max(R^2-x.^2,0).^2/24/delta;
elseif (strcmp(pot,'har')&& strcmp(initialcase,'lb'))
    mu=0.5*(3*gm*beta/2)^(2/3);
    rho=max((mu-gm^2*x.^2)/2,0)/beta;
else
    rho=exp(-gm*x.^2); dh=x(2)-x(1);
    lnorm=dh*sum(rho);
    rho=rho/lnorm;
end