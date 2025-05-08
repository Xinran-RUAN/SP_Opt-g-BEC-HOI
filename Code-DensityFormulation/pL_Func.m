function[pL_X] = pL_Func(X, L, data)
mass = 1;
domain_size = data.xmax - data.xmin;
G = GradObj_Func(X, data);
Xt = X - G / L;
pL_X = myPosConProj(Xt, mass, domain_size, 1e-12);
