% fraun_circ - Fraunholfer irradiance plot

L = 0.2;                %side length
M = 250;                % # samples
dx = L/M;               % sample interval
x = -L/2:dx:L/2-dx;
y = x;
[X,Y] = meshgrid(x,y);  %define the grid

w = 1e-3;               % x half-width
lamda = 0.633e-6;       % wavelength
z = 50;                 % propagation distance
k = 2*pi/lamda;         % wavenumber
lz = lamda*z;

% irradiance
I2 = (w^2/lz)^2.*(jinc(w/lz*sqrt(X.^2+Y.^2))).^2;

surf(x,y,nthroot(I2,3));
axis image;
