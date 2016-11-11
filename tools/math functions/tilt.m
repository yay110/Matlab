function u2 = tilt(u1,L,lambda,alpha,theta)

% tilt phasefront
% uniform sampling assumed
% u1 - input field
% L - side length
% lambda - wavelength
% alpha - tilt angle
% theta - rotaion angle (x axis)
% u2 - output field

[M,~] = size(u1);
dx = L/M;
k = 2*pi/lambda;

x = -L/2:dx:L/2-dx;
[X,Y]=meshgrid(x,x);

u2 = u1.*exp(1i*k*(X*cos(theta)+Y*sin(theta))...
    *tan(alpha));
end