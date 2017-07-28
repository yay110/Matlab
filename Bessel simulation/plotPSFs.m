for nn = 35

    fileName = ['psf',num2str(nn/10),'mm'];

    load(strcat(fileName,'.mat'));
    
    h = figure;
    normPsf = normr(Psf);
    Psf2F = Psf.^2;
    Psf3F = Psf.^3;
    displayLength2 = 30e-6;
    displayX = 300e-6;
    pixelNo = floor(displayLength2/dx2/2);
    displayCoord = M/2+1-pixelNo:M/2+1+pixelNo;
    subplot(3,1,1)
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,Psf(1:displayX*1e6,displayCoord)')
    axis image;
    xlabel('um');
    ylabel('um');
    title('Intensity of the PSF');
%     subplot(3,1,2)
%     imagesc(0:z:propagationDistance*1e6,x2(displayCoord)*1e6,normPsf(:,displayCoord)')
%     axis image;
%     xlabel('um');
%     ylabel('um');
%     title('Normalized Intensity of the PSF');
    subplot(3,1,2)
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,Psf2F(1:displayX*1e6,displayCoord)')
    axis image;
    xlabel('um');
    ylabel('um');
    title('2 photon Intensity of the PSF');
        subplot(3,1,3)
    imagesc(0:z:displayX*1e6,x2(displayCoord)*1e6,Psf3F(1:displayX*1e6,displayCoord)')
    axis image;
    xlabel('um');
    ylabel('um');
    title('3 photon Intensity of the PSF');
    
        saveas(h,strcat(fileName,'.svg'),'svg');
    saveas(h,strcat(fileName,'.png'),'png')
    
end