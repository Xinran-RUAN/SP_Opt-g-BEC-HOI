% test
x = (0:0.01:0.99)';
y = sin(2*pi*x);
N = length(x); 
Lx = 1;
mu = 2 * pi / Lx * (-N/2:N/2-1)';
Dy = D_x(y, mu);
plot(x, 2*pi*cos(2*pi*x), 'b'); hold on;
plot(x, Dy, 'r--'); hold off;
