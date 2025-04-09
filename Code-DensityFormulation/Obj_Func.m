function[E] = Obj_Func(X, data)
X = real(X); % check!
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;
N = data.Nx;
mu = 2 * pi / Lx * (-N/2:N/2-1)';
Dx = D_x(X, mu);

%% kinetic energy
E_kin = h / 8 * sum(abs(Dx).^2 ./ (X + vep)); 

%% potential energy
V = data.V;
E_pot = h * sum(V .* (sqrt(X.^2 + vep.^2) - vep));

%% beta - interaction
E_beta = 0.5 * beta * h * sum(X .^ 2);

%% delta - interaction
E_delta = 0.5 * delta * h * sum(abs(Dx).^2);

%%
E = E_kin + E_pot + E_beta + E_delta;
