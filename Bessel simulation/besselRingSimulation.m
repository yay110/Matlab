% lens optics Fourier transform
clear
tic;
% Define the sample grid
W1 = 0.36/5;       % sample region (m);
M = 10000;
dx1 = W1/M;
% x1 = -W1/2:dx1:W1/2-dx1;
x1 = -W1/2:dx1:W1/2;
[X1,Y1] = meshgrid(x1,x1);
[~, r] = cart2pol(X1,Y1);
% [theta,r] = cart2pol(x,y);

% define the optical system
lambda = 0.8e-6;
FL1 = 0.45;              % focal length of lens
FL2 = 0.009;
alpha = 2*pi/360;       % angle of the axicon is 1 degree
n = 1.4533;                 % refractive index of the axicon material
% NA = 0.5;

% define the incoming beam
w0 = 5e-3;      % radius of incoming beam
% g1 = exp(-(X1.^2+Y1.^2)/w0^2);                    %Amplitude of a Gaussian beam with radius w0.
g1 = exp(-r.^2/w0^2);
%
%  g1 = g1.*exp(-1i*alpha*sqrt(X1.^2+Y1.^2)*(n-1)/lambda*2*pi);
 g1 = g1.*exp(-1i*alpha*r*(n-1)/lambda*2*pi);
% g = exp(-(sqrt(X.^2+Y.^2)-ringRadius).^2/ringThickness^2);
Ig1 = abs(g1.^2);

% Display the original plane
displayLength1 = 20e-3;
pixelNo = displayLength1/dx1;
subplot(2,3,1)
imagesc(x1(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e3, x1(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e3,...
     Ig1((M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))));
axis image;
xlabel('mm');
ylabel('mm');
title(['pixel size is ', num2str(dx1*1e6),' um']);
 
 subplot(2,3,4);
plot(x1(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e3,Ig1(M/2+1,(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))));
xlabel('mm');


% ================ 1 st FFT=============
g1 = fftshift(g1);
g2 = fft2(g1);
g2  = fftshift(g2);
Ig2(:,:) = abs(g2).^2;

dx2 = 1/W1*FL1*lambda;
% x2 = -1/(2*dx1):1/W1:1/(2*dx1)-1/W1;
x2 = -1/(2*dx1):1/W1:1/(2*dx1);
x2 = x2 *FL1*lambda;
y2 = x2;
W2 = dx2 *M;


 displayLength2 = 10e-3;
 pixelNo = displayLength2/dx2;
 subplot(2,3,2);
imagesc(x2(M/2-floor(pixelNo/2)+1:M/2+1+floor(pixelNo/2))*1e3, x2(M/2+1-floor(pixelNo/2):M/2+1+floor(pixelNo/2))*1e3,...
     Ig2((M/2+1-floor(pixelNo/2):M/2+1+floor(pixelNo/2)),(M/2+1-floor(pixelNo/2):M/2+1+floor(pixelNo/2))));
 axis image;
 xlabel('mm');
 ylabel('mm');
 title(['pixel size is ', num2str(dx2*1e6),' um']);
 
  subplot(2,3,5);
plot(x2(M/2+1-floor(pixelNo/2):M/2+1+floor(pixelNo/2))*1e3,Ig2(M/2+1,(M/2+1-floor(pixelNo/2):M/2+1+floor(pixelNo/2))));
xlabel('mm');

%% =============== theory numbers comparison ===============

ringProfile = abs(g2(5001,5001:10001));
f = fit(x2(5001:10001)',ringProfile','gauss1');
plot(f,x2(5001:10001),ringProfile)
   thicknessRing = 2*f.c1;
   theoryThickness = 3.3*lambda*FL1/pi/w0;
   [~,radius] = max(ringProfile);
   radius = radius*dx2;
theoryRadius = (n-1)*alpha*FL1;

toc;