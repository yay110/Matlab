function createPressReleaseFigure()
    load 'z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-3\recording_lambda488nm_alpha3_beta100.mat';
    sample=max(restoredDataCube(1:400,1:275,1:500),[],3).';
    sample=sample./max(sample(:));
    dx=diff(xRange(1:2));
    dy=diff(yRange(1:2));
    xRange=dy*([1:275]-276/2);
    zRange=xRange;
    xRange=dx*([1:400]-400/2);
    
    objectiveNumericalAperture=0.4;
    pupilFunctor=@(U,V) exp(2i*pi*4*U.^3); %Bessel beam with 10% open fraction
    psfProj=calcVectorialPsf(yRange,zRange,xRange,500e-9,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,1.0,20,200e-3,[2]);
    psfProj=squeeze(psfProj);
    
    psfProj=psfProj./max(psfProj(:));
    img=0.1+repmat(sample,[1 1 3])+...
        cat(3,0*psfProj,psfProj,psfProj)+...
        3*cat(3,psfProj.*sample,0*psfProj.*sample,0*psfProj);
    %img=min(1,img);
    img=img./max(img(:));
    showImage(img,[],xRange,yRange);
    axis equal
end