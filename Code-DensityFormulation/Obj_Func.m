% X is real-valued
function[E] = Obj_Func(X, data)
beta = data.beta;
delta = data.delta;
vep =data.vep;
h = data.dx;
Lx = data.xmax - data.xmin;
N = data.Nx;

%% get diff_X = DX 
Lambda = 2 * pi * 1i / Lx * [0:(N/2-1), (-N/2:-1)];
Lambda = reshape(Lambda, size(X)); % Lambda与X形状保持一致
D_X = fourier_diff(X, Lambda);

%% kinetic energy
% E_kin = h / 8 * sum(abs(D_X).^2 ./ (X + vep)); 
E_kin = h / 8 * sum(abs(D_X).^2 ./ sqrt(X.^2 + vep^2)); 

%% potential energy
V = data.V;
E_pot = h * sum(V .* (sqrt(X.^2 + vep.^2) - vep));

%% beta - interaction
E_beta = 0.5 * beta * h * sum(X .^ 2);

%% delta - interaction
E_delta = 0.5 * delta * h * sum(abs(D_X).^2);

%%
E = E_kin + E_pot + E_beta + E_delta;
