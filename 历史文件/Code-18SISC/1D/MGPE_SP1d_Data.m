function data = MGPE_SP1d_Data(data)

if ~isfield(data,'beta'),    data.beta = 0; end
if ~isfield(data,'delta'),    data.delta = 0; end
if ~isfield(data,'init'),    data.init = 'Gaussian'; end
if ~isfield(data,'xmin'),    data.xmin = -16; end
if ~isfield(data,'xmax'),    data.xmax = 16; end
% if ~isfield(data,'Nx'),      data.Nx = 32*50+1; end
% if ~isfield(data,'vep'),      data.vep = 0.1; end

% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.x = (data.xmin:data.dx:data.xmax)';

% Potential function
x = data.x(2:end-1);
data.V = data.Potential(x);

