%
%
%
function createRayTracePlot(alpha,beta)
    close all;
    
    quickCalc=false;

    if nargin<1 || isempty(alpha)
        alpha=-3;
    end
    if nargin<2 || isempty(beta)
        beta=1; %0.3;
    end
    pupilAmplitude=@(U,V) (abs(U.^2+V.^2)>=(1-beta).^2).*(abs(U.^2+V.^2)<=1); %.*exp(-abs(U.^2+V.^2)/0.25);
    pupilOPDFunction=@(U,V) alpha*(U.^3 + V.^3);
    
    if beta<1
        outputFileName='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\setupIntro\BesselRayTrace';
    else
        if alpha==0
            outputFileName='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\setupIntro\GaussianRayTrace';
        else
            outputFileName='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\setupIntro\AiryRayTrace';
        end
    end
    
    absorberThickness=0.1;
    refractiveIndex=1.5;
    
    maskThickness=0.2;
    maskDrawnModulation=0.05;
    
    lensOPDFunction=@(U) -.5*abs(U).^2;
    lensPosition=0.4;
    lensThickness=0.5;
    focalDistance=1.5;
    
    nbRays=12;
    xLimits=[-2 5];
    
    uRange=[-1:.01:1];
    uRangeExt=[-1.5:.01:1.5];
    rayY=([1:nbRays]-1)/(nbRays/2); rayY=rayY-max(rayY)/2;
        
    fig=figure('Position',[50 50 800 300],'Color',[1 1 1],'NumberTitle','off','Name',outputFileName);
    ax=axes('Position',[0.1 0.1 0.8 0.8]);
    
    
    % Draw ray trace
    colormap(gray(256));
    f=fill([uRangeExt*0,uRangeExt*0-absorberThickness].',uRangeExt([1:end,end:-1:1]).',1+255*pupilAmplitude(uRangeExt([1:end,end:-1:1]),0).','FaceColor','interp','EdgeColor','none');
    set(f,'CDataMapping','direct'); hold on;
    
    for rayIdx=1:nbRays,
        xRange=[]; yRange=[];
        xRange(1:2)=[-0.5 -absorberThickness];
        yRange(1:2)=rayY(rayIdx);
        transmission=pupilAmplitude(yRange(1),0);
        if transmission>0
            raySlope=0;
            % Mask, back surface
            [xRange(3) yRange(3) raySlope]=findIntersection(@(U) maskThickness+maskDrawnModulation*pupilOPDFunction(U,0),1/refractiveIndex, xRange(2),yRange(2),raySlope);
            % Simulate perfect lens
            [xRange(4) yRange(4) raySlope]=paraxialLens(lensPosition+lensThickness,focalDistance,xRange(3),yRange(3),raySlope);
            xRange(5)=xRange(4); yRange(5)=yRange(4);
%             % Lens, front surface
%             [xRange(4) yRange(4) raySlope]=findIntersection(@(U) lensPosition-lensOPDFunction(U),refractiveIndex, xRange(3),yRange(3),raySlope);
%             % Lens, back surface
%             [xRange(5) yRange(5) raySlope]=findIntersection(@(U) lensPosition+lensThickness,1/refractiveIndex, xRange(4),yRange(4),raySlope);
            
            xRange(6)=xLimits(2);
            yRange(6)=yRange(5)+raySlope*diff(xRange(5:6));
        end
        p=plot(xRange,yRange,'LineWidth',2,'Color',[1 1 1]*.5);
        
    end
    
    fill([maskThickness+maskDrawnModulation*pupilOPDFunction(uRange,0), uRange*0],uRange([1:end,end:-1:1]),[0.8 0.8 1],'LineStyle','none');
    fill(lensPosition+[-lensOPDFunction(uRange), uRange*0+lensThickness],uRange([1:end,end:-1:1]),[0.8 0.8 1],'LineStyle','none');
    
    xlim(ax,xLimits); ylim(ax,[-2 2]);
    
    % Draw beam
    %
    wavelength=20e-3; %mm waves!
    refractiveIndex=1.0;
    objectiveNumericalAperture=1/focalDistance;
    resolution=.1*[[1 1]*wavelength/(2*objectiveNumericalAperture) wavelength*refractiveIndex/(objectiveNumericalAperture^2)];
    if quickCalc
        resolution(3)=resolution(3)*10;
    end
    xBeamRange=[-2:resolution(1):2];yBeamRange=[-5:resolution(2):5];zRange=[-focalDistance:resolution(3):(xLimits(2)-(lensPosition+lensThickness+focalDistance))];
    pupilFunctor=@(U,V) pupilAmplitude(U,V).*exp(2i*pi*pupilOPDFunction(U,V));%@(U,V) sqrt(U.^2+V.^2)>.9; %Bessel beam with 10% open fraction
    psfProj=calcVectorialPsf(xBeamRange,yBeamRange,zRange,wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,refractiveIndex,1,focalDistance*Inf,2);
    
    imAx=axes('Position',get(ax,'Position'));
    imAxPupil=axes('Position',get(ax,'Position'));
    psfProj=squeeze(psfProj)./max(psfProj(:));
    % psfProj=psfProj./repmat(max(psfProj),[size(psfProj,1) 1]); % normalize per slice
    %colorMapBeam=laserOnWhiteColorMap(1024,532e-9,1);
    colorMapBeam=interpolatedColorMap(1024,[1 1 1; 0 0 1; 0 1 .5; 1 0 0]);
    [XBeam,YBeam]=ndgrid(xBeamRange,yBeamRange);
    pupilAmplitudeProj=sum(pupilAmplitude(XBeam,YBeam),2); clear XBeam YBeam;
    pupilAmplitudeProj=pupilAmplitudeProj./max(pupilAmplitudeProj(:));
%     showImage(mapColor(repmat(pupilAmplitudeProj,[1 3]),colorMapBeam),[],-5+(0.5+[0:2])*(lensPosition+lensThickness--5)/3,xRange,imAxPupil); hold on;
    showImage(mapColor(psfProj(end:-1:1,:),colorMapBeam),[],pixelCentersToPixelEdges(lensPosition+lensThickness+focalDistance+zRange),pixelCentersToPixelEdges(xBeamRange),imAx); hold on;
    imH=get(imAx,'Children');
    %plainLaserColorImg=repmat(permute([0.25 1 0.25],[1 3 2]),[size(psfProj) 1]);
    plainLaserColorImg=mapColor(psfProj(end:-1:1,:),colorMapBeam);
    imwrite(plainLaserColorImg,[outputFileName,'.png'],'png','Alpha',1.5*psfProj(end:-1:1,:));
    % Swap color index for alpha value
    set(imH,'AlphaData',1.5*psfProj(end:-1:1,:));
    set(imH,'CData',plainLaserColorImg);
    xlim(imAx,xLimits);ylim(imAx,[-2 2]);
%     imagesc(zRange,xRange,squeeze(psfProj),'Parent',imAx);axis equal;xlabel('z');ylabel('x');
    axis(imAxPupil,'normal','off');
    axis(imAx,'normal','off');
    axis(ax,'normal','off');
    xlim(imAxPupil,xLimits);ylim(imAxPupil,[-2 2]);
    
    plot2svg([outputFileName,'.svg'],fig,'png');
end

function edgeRange=pixelCentersToPixelEdges(rng)
    dr=diff(rng(1:2));
    nbEls=numel(rng);
    edgeRange=rng(1)+(0.5+[1:nbEls]-1)*dr*(nbEls-1)/nbEls;
end

function [xPos yPos refractedRaySlope]=findIntersection(surfaceFunctor,refractiveIndexRatio, xPos,yPos,raySlope)
    % Calculate intersection
    precission=.01;
    maxIter=100;
    iter=0;
    xError=surfaceFunctor(yPos)-xPos;
    while iter<maxIter && abs(xError)>precission
        xPos=surfaceFunctor(yPos);
        yPos=yPos+xError*raySlope;
        xError=surfaceFunctor(yPos)-xPos; % update xError for next loop
        iter=iter+1;
    end
    if iter>=maxIter,
        logMessage('Ray trace did not converge in %d iterations',maxIter);
    end
    
    % Calculate refraction
    surfaceSlope=diff([0 0]+surfaceFunctor(yPos+[-.5 .5]/1000))*1000;
    surfaceSlopeWithRespectToRay=tan(atan(surfaceSlope)-atan(raySlope));
    refractedRaySlopeWithRespectToSurface=surfaceSlopeWithRespectToRay*refractiveIndexRatio;
    refractedRaySlope=tan(atan(surfaceSlope)-atan(refractedRaySlopeWithRespectToSurface));
end

function [xPos yPos newRaySlope]=paraxialLens(lensPosition,focalDistance,xPos,yPos,raySlope)
    % Trace until paraxial lens
    [xPos yPos]=findIntersection(@(U) lensPosition,1, xPos,yPos,raySlope);
    newRaySlope=raySlope-yPos/focalDistance;
end