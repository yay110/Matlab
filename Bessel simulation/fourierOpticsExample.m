% lens optics Fourier transform
clear
% Define the sample grid
W = 0.1;       % sample region (m);
M = 5000;
dx = W/M;
x = -W/2:dx:W/2-dx;
[X,Y] = meshgrid(x,x);
% [theta,r] = cart2pol(x,y);
thickness = [50,100,500,1000]*1e-6;


% scaling factor and scale;
lambda = 0.8e-6;
FL = 0.009;              % focal length of lens, for 20X, FL = 9mm.
NA = 0.5;
theta =asin(NA);
ringRadius = FL*tan(theta);

for n=1:4
% define the incoming beam
w0 = 5e-3;      % radius of incoming beam
% g = exp(-(X.^2+Y.^2)/w0^2);                    %Amplitude of a Gaussian beam with radius w0.
ringThickness = thickness(n);
g = exp(-(sqrt(X.^2+Y.^2)-ringRadius).^2/ringThickness^2);
Ig = abs(g.^2);
% g = g*exp(sqrt(-1)*((X./w0).^3+(Y./w0).^3));                      % Cubic phase

% displayLength = 20e-3;
% pixelNo = displayLength/dx;
% subplot(2,2,1)
% imagesc(x(M/2-pixelNo/2:M/2+pixelNo/2)*1e3, x(M/2-pixelNo/2:M/2+pixelNo/2)*1e3,...
%     Ig((M/2-pixelNo/2:M/2+pixelNo/2),(M/2-pixelNo/2:M/2+pixelNo/2)));
% axis image;
% xlabel('mm');
% ylabel('mm');
% title(['pixel size is ', num2str(dx*1e6),' um']);
% 
% subplot(2,2,2);
% plot(x(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e3,Ig(M/2+1,(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))));
% xlabel('mm');


% ================FFT=============
g0 = fftshift(g);
G0 = fft2(g0);
G  = fftshift(G0);
IG(:,:,n) = abs(G.^2);

dfx = 1/W*FL*lambda;
fx = -1/(2*dx):1/W:1/(2*dx)-1/W;
fx = fx *FL*lambda;
fy = fx;

% ==============display the graphs
 displayLength = 20e-6;
 pixelNo = displayLength/dfx;
% subplot(2,2,3);
% imagesc(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,fy(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,...
%     IG((M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))));
% axis image;
% xlabel('um');
% ylabel('um');
% title(['pixel size is ', num2str(dfx*1e9),' nm']);


% subplot(2,2,n);
% plot(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,IG(M/2+1,(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))));
% xlabel('um');
end

subplot(2,2,2);hold
for n = 1:4
plot(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,IG(M/2+1,(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),n)/max(max(IG(:,:,n))));
end
xlabel('um');legend('show');
title('Normalized intensity')

subplot(2,2,1);hold
for n = 1:4
plot(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,IG(M/2+1,(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),n));
end
xlabel('um');legend('show');
title('original intensity')

for n =1:4
    IGintegrated(:,n) = mean(squeeze(IG(:,:,n)));
end
    
subplot(2,2,3);hold
for n = 1:4
plot(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,IGintegrated((M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),n));
end
xlabel('um');egend('show');
title('original intensity - integrated to LS')

subplot(2,2,4);hold
for n = 1:4
plot(fx(M/2-floor(pixelNo/2):M/2+floor(pixelNo/2))*1e6,IGintegrated((M/2-floor(pixelNo/2):M/2+floor(pixelNo/2)),n)/max(IGintegrated(:,n)));
end
xlabel('um');legend('show');
title('normalized intensity - integrated to LS')
