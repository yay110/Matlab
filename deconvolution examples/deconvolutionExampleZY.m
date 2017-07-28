% function deconvolutionExample2D()


I = imread('cameraman.tif');
imgSize = size(I);
%     imgSize = [100, 100];
xRange = [1:imgSize(1)]-1-floor(imgSize(1)/2);
yRange = [1:imgSize(2)]-1-floor(imgSize(2)/2);
%   create the grid
[X,Y] = ndgrid(xRange,yRange);

obj = double(I)./max(max(double(I)));
%     obj = @(x,y) ((x==10&y==10)|(x==-25&y==-25));
%     obj = @(x,y) (x==50&&y==60)

Psf = exp(-X.^2./(2*5^2)).*exp(-Y.^2./(2*5^2));


OTF=fft2(ifftshift(Psf));
OBJ=fft2(obj);
%     OBJ=fft2(obj(X,Y));
IMG=OTF.*OBJ;
Img=ifft2(IMG);
Img=Img+0.01*randn(imgSize)*max(max(Img));

subplot(2,3,1);imagesc(obj);colorbar;
subplot(2,3,2);imagesc(Psf);
subplot(2,3,3);imagesc(Img);colorbar;
% clear obj OBJ IMG;

dx=diff(xRange(1:2));
dy=diff(yRange(1:2));
dfx=1./(dx.*imgSize(1));
dfy=1./(dy.*imgSize(2));
fxRange = ([1:imgSize(1)]-1-floor(imgSize(1)/2))*dfx;
fyRange = ([1:imgSize(2)]-1-floor(imgSize(2)/2))*dfy;
[FX,FY] = ndgrid(fxRange,fyRange);

%  freqency domain, normalized by the diffraction limit which is assumed to
%  be 1/3 of the highest frequency this system can acquire. 
fRel = ifftshift(sqrt(FX.^2+FY.^2)./(max(FX(:))/3));

%   This is assuming that the noise increases with frequency. 
NSR = fRel./3;

%     Hw = conj(OTF)./(conj(OTF).*OTF+NSR.^2);
Hw = conj(OTF)./(abs(OTF).^2+NSR.^2);


IMG = fft2(Img);
Img_decon = IMG.* Hw;
Img = ifft2(Img_decon);

subplot(2,3,4);imagesc(abs(fftshift(IMG)));
subplot(2,3,5);imagesc(max(0,Img));colorbar;









% end