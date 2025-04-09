function[G] = GradObj_Func(X, data)
X = real(X); % check!
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;
N = data.Nx;
mu = 2 * pi / Lx * (-N/2:N/2-1)';
Dx = D_x(X, mu);
cDx = conjD_x(X, mu);
%% kinetic energy
G_kin = h * (tranD_x(cDx./(X+vep), mu) + conjtranD_x(Dx./(X+vep), mu) - abs(Dx).^2./(X+vep).^2);

%% potential energy
V = data.V;
G_pot = h * V .* X ./ sqrt(X.^2 + vep.^2);

%% beta - interaction
G_beta = h * beta * X;

%% delta - interaction
G_delta = 0.5 * delta * h * (tranD_x(cDx, mu) + conjtranD_x(Dx, mu));


% figure
% plot(G_kin, 'b'); hold on;
% plot(G_pot, 'r'); hold on;
% plot(G_beta, 'k'); hold on;
% plot(G_delta, 'g'); hold off;
%%
G = G_kin + G_pot + G_beta + G_delta;