function w = Prox(v, t, opt)
b=opt.b; data=opt.data; N=data.Nx-2; 
v=v-min(v)+2*b/N;
w = ProjectOntoL1Ball(v, b);
