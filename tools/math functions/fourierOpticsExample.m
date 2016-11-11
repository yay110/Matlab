% lens optics Fourier transform
clear
% Define the sample grid
W = 0.01;       % sample region (m);
M = 200;
dx = W/M;
x = -W/2:dx:W/2-dx;
[X,Y] = meshgrid(x,x);
% [theta,r] = cart2pol(x,y);

% define the incoming beam
w0 = 0.001;      % radius of incoming beam
e = exp(-(X.^2+Y.^2)/w0^2);                    %Amplitude of a Gaussian beam with radius w0.

subplot(2,2,1)
imagesc(x*1e3,x*1e3,abs(e))
axis image;
xlabel('mm');
ylabel('mm');

g0 = fftshift(e);
G0 = fft2(g0);
G  = fftshift(G0);

% scaling factor
lambda = 0.5e-6;
F = 0.1;

fx2 = linspace(-lambda*F*M/(2*W), lambda*F*M/(2*W), M);


fx = -1/(2*dx):1/W:1/(2*dx)-1/W;
fx = fx *F*lambda;


fy = fx;

subplot(2,2,3);
imagesc(fx*1e6,fy*1e6,abs(G));
axis image;
xlabel('um');
ylabel('um');

subplot(2,2,4);
plot(fx*1e6,abs(G(M/2+1,:)));
xlabel('um');