function data = MGPE_FD2d_Data(data)

if ~isfield(data,'beta'),    data.beta = 0; end
if ~isfield(data,'delta'),    data.delta = 0; end
if ~isfield(data,'init'),    data.init = 'Gaussian'; end
if ~isfield(data,'xmin'),    data.xmin = -16; end
if ~isfield(data,'xmax'),    data.xmax = 16; end
if ~isfield(data,'ymin'),    data.ymin = -16; end
if ~isfield(data,'ymax'),    data.ymax = 16; end


% Computational domain
data.dx = (data.xmax-data.xmin)/(data.Nx-1);
data.dy = (data.ymax-data.ymin)/(data.Ny-1);
data.x = (data.xmin:data.dx:data.xmax)';
data.y = (data.ymin:data.dy:data.ymax)';

% Potential function
x = data.x(2:end-1); y = data.y(2:end-1);
[X,Y]=ndgrid(x,y);
data.V = data.Potential(X,Y);

