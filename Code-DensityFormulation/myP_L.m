function[pL_X] = myP_L(X, L, data)
mass = 1;
domain_size = data.xmax - data.xmin;
[~, G] = myDiscreteEnergy(X, data);
Xt = X - G / L;
pL_X = myPosConProj(Xt, mass, domain_size, 1e-12);
