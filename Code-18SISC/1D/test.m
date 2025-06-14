h=0.5; % change h in myfun!
xx_full=(-16:h:16)'; xx=xx_full(2:end-1);
X0 = exp(-xx.^2); X0 = X0/norm(X0,1)*1/h;

A=[];b=[];lb=zeros(size(X0));ub=[];nonlcon = [];
Aeq=ones(1,length(X0)); beq=1/h;
options.MaxFunEvals=100000; 

x = fmincon(@myfun,X0,A,b,Aeq,beq,lb,ub,nonlcon,options);
myfun(x)