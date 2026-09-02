function[fx]=myfun(x)
vep=0.1; h=0.5; beta=0;delta=0;
xx_full=(-16:h:16)'; xx=xx_full(2:end-1);
V=0*xx.^2/2;
fx=1/2/h*((sqrt(x(1)+vep)-sqrt(vep))^2+sum((diff(sqrt(x+vep))).^2)+(sqrt(vep)-sqrt(x(end)+vep))^2)...
    +h*sum(V.*x)+h*beta/2*sum(x.^2)+delta/2/h*(sum(diff(x).^2)+(x(1)^2+x(end)^2));
