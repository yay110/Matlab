% define the optical system
lambda = 0.8e-6;
FL1 = 0.3;               % focal length of the first lens, after the axicon
FL2 = 0.018;            % focal length of the objective lens
alpha = 2*pi/360;       % angle of the axicon is 1 degree
n = 1.4533;             % refractive index of the axicon material
% NA = 0.5;

% =========================== define grid  ========================
W1 = 0.1;               % sample region (m);
M = 10000;
dx1 = W1/M;             % pixel size
x1 = -W1/2:dx1:W1/2;    % define coordinate
[X,Y] = meshgrid(x1,x1);
[~, r1] = cart2pol(X,Y);

w0=5e-3;                % radius of the incoming Bessel beam

for nn =5
    
    w0 = nn/2*1e-3;
    % ============= calculate the parameter of ring at back apeture===========
    % based on paper Pierre-Andre Belanger 1978.
    theoryThickness = 3.3*lambda*FL1/pi/w0;     % Thickness of the ring
    theoryRadius = (n-1)*alpha*FL1;             % diameter of the ring
    n1=1;
    trueNA = sin(atan(theoryRadius/FL2));       % Resulting NA with the diametre of ring
    beta = theoryThickness/theoryRadius;        % Size of the ring
    theoryFOV = lambda/n1 / (2*(1-sqrt(1-(trueNA/n1)^2))*beta); % FOV based on Airy paper
    
    % ======================== beam at back aperture ====================
    g = exp(-(r1-theoryRadius).^2/(theoryThickness/2)^2);
    Ig = abs(g).^2;
    
    % ======================= display the original plane ======================
    displayLength = 20e-3;
    pixelNo = displayLength/dx1;
    figure
    subplot(2,2,1);
    displayCoord = M/2+1-pixelNo/2:M/2+1+pixelNo/2;
    imagesc(x1(displayCoord)*1e3, x1(displayCoord)*1e3,Ig(displayCoord,displayCoord));
    axis image;
    xlabel('mm');
    ylabel('mm');
    title('Intensity of the ring on back aperture.');
    subplot(2,2,3);
    plot(x1(displayCoord)*1e3,Ig(M/2+1,displayCoord));
    xlabel('mm');
    title('Profile of the ring on back aperture.')
    
    % ================ FFT to focal plane    =============
    g = fftshift(g);
    g2 = fft2(g);
    g2  = fftshift(g2);
    Ig2 = abs(g2).^2;
    
    dx2 = 1/W1*FL2*lambda;
    x2 = -1/(2*dx1):1/W1:1/(2*dx1);
    x2 = x2 *FL2*lambda;
    y2 = x2;
    W2 = dx2 *M;
    
    
    % ============== display the focal plane Bessel beam ==================
    displayLength2 = 20e-6;
    pixelNo = floor(displayLength2/dx2/2);
    displayCoord = M/2+1-pixelNo:M/2+1+pixelNo;
    subplot(2,2,2);
    imagesc(x2(displayCoord)*1e6, x2(displayCoord)*1e6,...
        Ig2(displayCoord,displayCoord));
    axis image;
    xlabel('um');
    ylabel('um');
    title('Intensity of the Bessel beam on the focal plane');
    subplot(2,2,4);
    plot(x2(displayCoord)*1e6,Ig2(M/2+1,displayCoord));
    xlabel('um');
    title('Profile of the  Bessel beam on focal plane');
    
    
    % =========== calculating effeciency ========================
    centreIntensity = Ig2(M/2+1,:);
    psf = centreIntensity(1,M/2+1:end);
    [~,localM] = findpeaks(-1*psf);
    [rr,cc] = meshgrid(1:M+1,1:M+1);
    C = sqrt((rr-M/2-1).^2+(cc-M/2-1).^2)<=localM(1);
    centreCore = Ig2.*C;
    effi(nn) = sum(centreCore(:))/sum(Ig2(:));
    %Circle with radius centered at 40,40 in an 80x80image
%     effiLine = sum(centreIntensity(M/2-localM(1)+2:M/2+localM(1)))/sum(centreIntensity(:));
end