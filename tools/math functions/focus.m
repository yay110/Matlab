function u2 = focus(u1,L,lambda,zf)

% converging or diverging phase-front
% uniform sampling assumed
% u1 - input field
% L - side length
% lambda - wavelength
% zf - focal distance( + converging, - diverging)
% u2 - output field

[M,~] = size(u1);       % get input field array size
dx = L/M;               % sample interval
k = 2*pi/lambda;

x = -L/2:dx:L/2-dx;
[X,Y]=meshgrid(x,x);

u2 = u1.*exp(-1i*k/(2*zf)*(X.^2+Y.^2));     % apply focus
end