%
%  Program Axicon - 8/4/2013
%  Solves field propagation past an axicon
%
%function Axicon;
%
%  Read in input data
%
clear; 

w0 = 2e3;              %unit micron, radius of incoming beam
lambda = 0.8;           %unit micron
k = 2*pi/lambda;        
axang = 178;            %axcon angle
n = 1.4533;                %refractive index of axicon
%
beamAngle = (n-1)*(180-axang)*pi/360;

%% ================ Grid Parameters =================
% xymax   = 10*1e3;          %simulation range for XY is 10mm cross
% nxy     = 512;
% zmax    = 100*1e3;         %simulation range for Z is 100mm long
% nz      = 512;
kr = k*(n-1)*(180-axang)*pi/360;       %I believe this is not used

xymax   = 10000;
nxy     = 1000;
zmax    = w0/(kr/k);         % length of Bessel beam zmax
L       = 3*zmax;
nz      = 100;              %Number of points in length of L
R       = 100e3;            %Define the focal length of the spherical lens
%
dz = L/nz;
dz = 0.5e3;
xmax = xymax;
dx = xmax/nxy;
kmax = 2*pi/dx;
dk = kmax/nxy;
nmid = floor(nxy/2);

v = [0:nxy-1];
[x,y] = meshgrid(v,v);
x = x*dx - xmax/2;
y = y*dx - xmax/2;
xyar = v*dx - xmax/2;
p = find(v > nmid);
v(p) = nxy-v(p);
v = v*dk;
[k2,k1] = meshgrid(v,v);
[theta,r] = cart2pol(x,y);

% Parameter for diffraction simulation
arg   = -dz*(k1.^2+k2.^2)/2/k;
freq  = exp(sqrt(-1)*arg);

%
%%================ Set up input beam ===========================
%
e = exp(-(r-2000).^2/50^2);

phaseAxicon = -r*(180-axang)*pi/360*(n-1)/lambda; 
e = exp(-r.^2/w0^2);                    %Amplitude of a Gaussian beam with radius w0.
e = e.*exp(1i*phaseAxicon*2*pi);        %adding the phase of axicon
% e = e.*exp(-1i*k*r.^2/2/R);             %adding the phase of a shperical lens

f = fft2(e);
f = fftshift(f);

e = e.*exp(-1i*k*r.^2/2/R);             %adding the phase of a shperical lens

z = 0;
zval(1)    = 0;
ii=1;
iprof(ii,:) = abs(e(:,nmid+1)).^2;

for i = 1:3*nz
    ii=ii+1;
    e = fft2(ifft2(e).*freq);
    imax = max(abs(e(:,nmid+1)).^2);
    iprof(ii,:) = abs(e(:,nmid+1)).^2/imax;
%     iprof(ii,:) = abs(e(:,nmid+1));
    zval(ii) = zval(ii-1) + dz;
end

subplot(1,2,1)
imagesc(xyar./1e3,xyar./1e3,abs(f));axis image;
subplot(1,2,2)
imagesc(zval./1e3,xyar./1e3,iprof');
set(gca,'FontSize',8);
xlabel('z (mm)');
ylabel('x (mm)');
title('Normalized intensity |E(x,0,z)| ^2');

% e = e.*exp(1i*phaseAxicon*2*pi);        %adding the phase of axicon

% for i = 1:20
%     ii=ii+1;
%     e = fft2(ifft2(e).*freq);
%     imax = max(abs(e(:,nmid+1)).^2);
%     iprof(ii,:) = abs(e(:,nmid+1)).^2/imax;
% %     iprof(ii,:) = abs(e(:,nmid+1));
%     zval(ii) = zval(ii-1) + dz;
% end
% 
% e = e.*exp(-1i*k*r.^2/2/R);
% 
% for i = 1:5*nz
%     ii=ii+1;
%     e = fft2(ifft2(e).*freq);
%     imax = max(abs(e(:,nmid+1)).^2);
%     iprof(ii,:) = abs(e(:,nmid+1)).^2/imax;
% %     iprof(ii,:) = abs(e(:,nmid+1));
%     zval(ii) = zval(ii-1) + dz;
% end
%

figure;hold on;
subplot(2,5,[1 2 3 4 6 7 8 9]);
 imagesc(zval./1e3,xyar./1e3,iprof');
set(gca,'FontSize',8);
xlabel('z (mm)');
ylabel('x (mm)');
title('Normalized intensity |E(x,0,z)| ^2');
subplot(2,5,5);
imagesc(xyar((256-140):(256+140)),xyar((256-140):(256+140)),abs(e((256-140):(256+140),(256-140):(256+140))).^2);axis image;axis off;
xlabel('x (\mum)');
ylabel('y (\mum)');
subplot(2,5,10);
imagesc(xyar((256-140):(256+140)),xyar((256-140):(256+140)),angle(e((256-140):(256+140),(256-140):(256+140))));axis image;axis off;
xlabel('x (\mum)');
ylabel('y (\mum)');
%
% figure(2)
% %plot(xyar,abs(e(nmid+1,:)).^2,xyar,iprof(end,:));
% plot(xyar,iprof(1,:),xyar,iprof(end,:));
% set(gca,'FontSize',15);
% xlabel('x (microns)');
% ylabel('|E(x,0,z=L)|^2');
% title('Final intensity profile');
% %
% figure(3)
% subplot(2,2,1)
% imagesc(xyar,xyar,abs(e).^2);
% colorbar
% set(gca,'FontSize',15);
% xlabel('x (microns)');
% ylabel('y (microns)');
% title('Final intensity profile');
% %
% subplot(2,2,2)
% imagesc(xyar,xyar,angle(e));
% colorbar
% set(gca,'FontSize',15);
% xlabel('x (microns)');
% ylabel('y (microns)');
% title('Final Phase profile');