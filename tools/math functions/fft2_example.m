%fft2_example - 2D FFT example

% Define the signal
wx = 0.1;
wy = 0.05;
L  = 2;
M  = 200;
dx = L/M;
x  = -L/2:dx:L/2;
y  = x;
[X,Y] = meshgrid(x,y);
g = rect(X/(2*wx)).*rect(Y/(2*wy));
subplot(2,2,1);
mesh(x,y,g);

% Perfrom a 2D fft of the rectangle function
g0 = fftshift(g);
G0 = fft2(g0)*dx^2;
G  = fftshift(G0);
fx = -1/(2*dx):1/L:1/(2*dx);
fy = fx;

subplot(2,2,3);
surf(fx,fy,abs(G));
camlight left; lighting phong;
shading interp;
ylabel('fy (cyc/m)');xlabel('fx (cyc/m)');

subplot(2,2,4);
plot(fx,abs(G(M/2+1,:)));
title('Magnitude');
xlabel('fx (cyc/m)');