clear;
% define the optical system
lambda = 0.8e-6;
FL1 = 0.45;             % focal length of the first lens, after the axicon
FL2 = 0.018;            % focal length of the objective lens
alpha = 2*pi/360;       % angle of the axicon is 1 degree
n = 1.4533;             % refractive index of the axicon material
% NA = 0.5;
w0=2.5e-3;                % radius of the incoming Bessel beam

for nn = 5
    
    w0 = nn/2*1e-3;
    
    % ============= calculate the parameter of ring at back apeture===========
    % based on paper Pierre-Andre Belanger 1978.
    theoryThickness = 3.3*lambda*FL1/pi/w0;     % Thickness of the ring
    theoryRadius = (n-1)*alpha*FL1;             % diameter of the ring
    n1=1;
    trueNA = sin(atan(theoryRadius/FL2));       % Resulting NA with the diametre of ring
    beta = theoryThickness/theoryRadius;        % Size of the ring
    theoryFOV = lambda/n1 / (2*(1-sqrt(1-(trueNA/n1)^2))*beta); % FOV based on Airy paper
    
    % =========================== define grid  ========================
    W1 = 0.1;               % sample region (m);
    M = 5000;
    dx1 = W1/M;             % pixel size
    x1 = -W1/2:dx1:W1/2;    % define coordinate
    [X,Y] = meshgrid(x1,x1);
    [~, r1] = cart2pol(X,Y);
    
    % ======================== beam at back aperture ====================
    g = exp(-(r1-theoryRadius).^2/(theoryThickness/2)^2);
    Ig = abs(g).^2;
    
    % ======================= display the original plane ======================
    %     displayLength = 20e-3;
    %     pixelNo = displayLength/dx1;
    %     figure
    %     subplot(2,2,1);
    %     displayCoord = M/2+1-pixelNo/2:M/2+1+pixelNo/2;
    %     imagesc(x1(displayCoord)*1e3, x1(displayCoord)*1e3,Ig(displayCoord,displayCoord));
    %     axis image;
    %     xlabel('mm');
    %     ylabel('mm');
    %     title('Intensity of the ring on back aperture.');
    %     subplot(2,2,3);
    %     plot(x1(displayCoord)*1e3,Ig(M/2+1,displayCoord));
    %     xlabel('mm');
    %     title('Profile of the ring on back aperture.')
    
    
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
    %     displayLength2 = 20e-6;
    %     pixelNo = floor(displayLength2/dx2/2);
    %     displayCoord = M/2+1-pixelNo:M/2+1+pixelNo;
    %     subplot(2,2,2);
    %     imagesc(x2(displayCoord)*1e6, x2(displayCoord)*1e6,...
    %         Ig2(displayCoord,displayCoord));
    %     axis image;
    %     xlabel('um');
    %     ylabel('um');
    %     title('Intensity of the Bessel beam on the focal plane');
    %     subplot(2,2,4);
    %     plot(x2(displayCoord)*1e6,Ig2(M/2+1,displayCoord));
    %     xlabel('um');
    %     title('Profile of the  Bessel beam on focal plane');
    
    %================ propagation simulation ==========================
    
    fx2 = -1/(2*dx2):1/W2:1/(2*dx2);    % frequence coords
    [FX2,FY2] = meshgrid(fx2,fx2);
    [~,r2] = cart2pol(FX2,FY2);
    
    z = 10e-6;                           % propogation step
    propagationDistance = 500e-6;       % simulation distance
    
    H = exp(-1i*pi*lambda*z*r2.^2);     % transfer function
    H = fftshift(H);
    G2 = fft2(fftshift(g2));
    centreLine = g2(M/2+1,:);
    Psf(1,:) = abs(centreLine).^2;
    
    integratedPsf(1,:) = mean(abs(g2).^2);
    
    
    % ==================== propagation ================================
    for n=1:propagationDistance/z
        n
        U2 = H.*G2;
        u2 = ifftshift(ifft2(U2));
        centreLine = u2(M/2+1,:);               % Centre line of the field
        Psf(n+1,:) = abs(centreLine).^2;        % Intensity of the centre Line to PSF
        integratedPsf(n+1,:) = mean(abs(u2).^2);
        G2 = U2;
    end
    
    h = figure;
    %     normPsf = normr(Psf);
    Psf2F = integratedPsf.^2;
    Psf3F = integratedPsf.^3;
    displayLength2 = 50e-6;
    pixelNo = floor(displayLength2/dx2/2);
    displayCoord = M/2+1-pixelNo:M/2+1+pixelNo;
    displayX = 300e-6;
    subplot(3,3,[1,2])
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,integratedPsf(1:(displayX/z),displayCoord)')
    axis image;
    xlabel('\mum');
    ylabel('\mum');
    title('Single photon Intensity of the PSF - integrated');
    subplot(3,3,3)
    plot(x2(displayCoord)*1e6,integratedPsf(1,displayCoord)/max(integratedPsf(:)));
    view(90,90) % swap the x and y axis
    ylabel('Normalized intensity');
    xlabel('\mum');
    title('Intensity profile at 0 (focus)');
    
    subplot(3,3,[4,5])
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,Psf2F(1:(displayX/z),displayCoord)')
    axis image;
    xlabel('\mum');
    ylabel('\mum');
    title('2 photon Intensity of the PSF - integrated');
    subplot(3,3,6)
    plot(x2(displayCoord)*1e6,Psf2F(1,displayCoord)/max(Psf2F(:)));
    view(90,90) % swap the x and y axis
    ylabel('Normalized intensity');
    xlabel('\mum');
    title('Intensity profile at 0 (focus)');
    subplot(3,3,[7,8])
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,Psf3F(1:(displayX/z),displayCoord)')
    axis image;
    xlabel('\mum');
    ylabel('\mum');
    title('3 photon Intensity of the PSF - integrated');
    subplot(3,3,9)
    plot(x2(displayCoord)*1e6,Psf3F(1,displayCoord)/max(Psf3F(:)));
%     view(90,90) % swap the x and y axis
    ylabel('Normalized intensity');
    xlabel('\mum');
    title('Intensity profile at 0 (focus)');
    
    
    
    %     fileName = ['psf',num2str(w0*1e3),'mm'];
    %     save(strcat(fileName,'.mat'),'Psf','dx2','M','z','propagationDistance','x2');
    %     saveas(h,strcat(fileName,'.svg'),'svg');
    %     saveas(h,strcat(fileName,'.png'),'png')
end

