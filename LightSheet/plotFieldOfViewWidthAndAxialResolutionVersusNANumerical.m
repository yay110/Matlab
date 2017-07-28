% [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNANumerical()
%
% Calculates the field of view and axial resolution of a light sheet microscope as a function of the illumination NA.
%
% All SI units.
%
function [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNANumerical()
    close all;
    
    config=struct();
    config.excitation=struct();
    config.excitation.objective=struct();
    config.excitation.objective.refractiveIndex=1.00;
    config.excitation.objective.numericalAperture=0.42;
    config.excitation.objective.magnification=20;
    config.excitation.objective.tubeLength=0.200;
    config.excitation.wavelength=532e-9;
	config.excitation.fractionOfNumericalApertureUsed=1;
    config.sample.refractiveIndex=1.40;
    config.detection=struct();
    config.detection.objective=struct();
    config.detection.objective.refractiveIndex=1.00;
    config.detection.objective.numericalAperture=0.40;
    config.detection.objective.magnification=20;
    config.detection.objective.tubeLength=0.160;
    config.detection.tubeLength=0.176;
    config.detection.wavelength=600e-9; %6.12e-7;
	config.detection.fractionOfNumericalApertureUsed=1;
    config.detector=struct();
    config.detector.pixelSize=[1 1]*7.4e-6;   
    
    stepSizeNA=0.001;
    numericalApertures=stepSizeNA:stepSizeNA:(config.sample.refractiveIndex);
    
    alpha=7;
    betas=[0.10 0.05];
    
    calcToShow=2; % 1: theoretical, 2: numerical
    
    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    depthOfFocus=config.sample.refractiveIndex*(config.detection.wavelength/(config.detection.objective.numericalAperture^2)+...
        config.detector.pixelSize(1)/(realMagnification*config.detection.objective.numericalAperture));
    spotLength=4*config.sample.refractiveIndex*config.detection.wavelength/(config.detection.objective.numericalAperture^2);
    
    fontSize=24;
    fontName='Arial';
    lineColors=[.33 .33 .33; 0.2 0.2 0.8; 0.8 0 0; 0 .5 0];
    %lineStyles={':','-','--','-'};
    lineStyles={'-','-','-','-'};
    lineWidths=[2 3 2 3];
    
    % Calculate theoretical values
    FOVGaussian=4*config.excitation.wavelength*config.sample.refractiveIndex*numericalApertures.^-2;
    axialResolutionGaussian(1,:)=config.excitation.wavelength./(2*numericalApertures*0.88);
    FOVBessel10=(config.excitation.wavelength/config.sample.refractiveIndex)./(2*(1-sqrt(1-(numericalApertures/config.sample.refractiveIndex).^2))*betas(1));
    axialResolutionBessel10(1,:)=config.excitation.wavelength*pi*.05./(numericalApertures*betas(1));
    FOVBessel5=(config.excitation.wavelength/config.sample.refractiveIndex)./(2*(1-sqrt(1-(numericalApertures/config.sample.refractiveIndex).^2))*betas(2));
    axialResolutionBessel5(1,:)=config.excitation.wavelength*pi*.05./(numericalApertures*betas(2));
    FOVAiry=(6*alpha*config.excitation.wavelength/config.sample.refractiveIndex)./(1-sqrt(1-(numericalApertures/config.sample.refractiveIndex).^2));
    maxSpFreqAiry=min(0.88,1/(0.05^2*48*alpha)).*numericalApertures/(config.excitation.wavelength/2);
    axialResolutionAiry(1,:)=1./maxSpFreqAiry;
    
    % Calculate the numerical values
    for numericalApertureIdx=1:length(numericalApertures)
        numericalAperture=numericalApertures(numericalApertureIdx);
        logMessage('Calculating axial resolution for NA=%0.3f',numericalAperture);
        axialResolutionGaussian(2,numericalApertureIdx)=calcAxialResolution(config,numericalAperture,0,1);
        axialResolutionBessel10(2,numericalApertureIdx)=calcAxialResolution(config,numericalAperture,0,.10);
        axialResolutionBessel5(2,numericalApertureIdx)=calcAxialResolution(config,numericalAperture,0,.05);
        axialResolutionAiry(2,numericalApertureIdx)=calcAxialResolution(config,numericalAperture,7,1);
    end
    
    save('FieldOfViewWidthAxialResolutionNA4.mat','numericalApertures','axialResolutionGaussian','axialResolutionBessel10','axialResolutionBessel5','axialResolutionAiry','FOVGaussian','FOVBessel10','FOVBessel5','FOVAiry','depthOfFocus','spotLength');
    
    if (nargout==0)
        % Display results
        figs(1)=figure('Position',[100 100 1024 768]);
        axs(1)=subplot(2,2,1,'Parent',figs(1));
        axs(2)=subplot(2,2,2,'Parent',figs(1));
        axs(3)=subplot(2,2,[3 4],'Parent',figs(1));
        labels=[];
        semilogy(numericalApertures,axialResolutionGaussian(calcToShow,:)*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionBessel10(calcToShow,:)*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionBessel5(calcToShow,:)*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionAiry(calcToShow,:)*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionGaussian(calcToShow,:)*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(1)); hold(axs(1),'on');
%         semilogy(numericalApertures,spotLengths*1e6,':','LineWidth',2,'Color',[1 1 1]*.33,'Parent',axs(1)); hold(axs(1),'on');
        xlim(axs(1),[0 numericalApertures(end)]); ylim(axs(1),[.1 500]);
        set(axs(1),'XTick',[0 0.33 0.67 1 1.33]); set(axs(1),'YTick',10.^[0:1:100]);
        labels(end+1)=xlabel(axs(1),'NA');
        labels(end+1)=ylabel(axs(1),'Axial Resolution [\mum]');
        legend(axs(1),{'Gaussian','Bessel10','Bessel5','Airy'});
        semilogy(numericalApertures,FOVGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(2)); hold(axs(2),'on');
        xlim(axs(2),[0 numericalApertures(end)]); ylim(axs(2),[.1 500]);
        set(axs(2),'XTick',[0 0.33 0.67 1 1.33]); set(axs(2),'YTick',10.^[0:1:100]);
        labels(end+1)=xlabel(axs(2),'NA');
        labels(end+1)=ylabel(axs(2),'FOV Width [\mum]');
        legend(axs(2),{'Gaussian','Bessel10','Bessel5','Airy'});
%         loglog(axialResolutionGaussian(1,:)*1e6,FOVGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionBessel10(1,:)*1e6,FOVBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',[0.8 0.2 0.2],'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionBessel5(1,:)*1e6,FOVBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionAiry(1,:)*1e6,FOVAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVGaussian*1e6,axialResolutionGaussian(calcToShow,:)*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVBessel10*1e6,axialResolutionBessel10(calcToShow,:)*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVBessel5*1e6,axialResolutionBessel5(calcToShow,:)*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVAiry*1e6,axialResolutionAiry(calcToShow,:)*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(3)); hold(axs(3),'on');
        xlim(axs(3),[0 1000]); ylim(axs(3),[0 10]);
        set(axs(3),'XTick',[0:100:10000]); set(axs(3),'YTick',[0:2:100]);
        labels(end+1)=xlabel(axs(3),'FOV Width [\mum]');
        labels(end+1)=ylabel(axs(3),'Axial Resolution [\mum]');
        legend(axs(3),{'Gaussian','Bessel10','Bessel5','Airy'},'Parent',figs(1));
        set(axs,'LineWidth',3,'FontSize',fontSize,'FontSize',fontSize,'FontWeight','bold');
        set(labels,'FontSize',fontSize,'FontWeight','bold');
                
        plot2svg('C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\resolutionVsFOV\resolutionVsFOVNumerical_raw.svg',figs(1));
        
        clear numericalApertures;
    end
end
function axialResolution=calcAxialResolution(config,numericalAperture,alpha,beta)
    config.excitation.objective.numericalAperture=numericalAperture;

    mtfThreshold=0.05;
    
    % Slightly oversample
    zStepSize=0.9*0.5*config.excitation.wavelength/(2*config.excitation.objective.numericalAperture);
    nbSamplesZ=1024;
    
    xRange=0;
    yRange=0;
    zRange=([1:nbSamplesZ]-floor(nbSamplesZ/2)-1)*zStepSize;

    logMessage('Calculating light sheet...');
    tilt=-alpha*3/5;
    lightSheetPsf=calcLightSheetPsf(single(xRange),single(yRange),single(zRange),tilt,config.excitation,alpha,beta,config.sample.refractiveIndex);
    lightSheetPsf=lightSheetPsf(:).';
    
    % Calc PSF projection and OTF slice
    otfStep=(1./diff(zRange(1:2)))/length(zRange);
    ZOtf=([0:length(zRange)-1]-floor(length(zRange)/2))*otfStep;
    otf=fftshift(fft(ifftshift(lightSheetPsf)));
    mtf=abs(otf); mtf=mtf./mtf(1+floor(end/2));
    
    firstIndexBelow5pct=find(ZOtf>=0 & mtf<mtfThreshold,1,'first');
    if (~isempty(firstIndexBelow5pct))
        % Linearly interpolate
        spFreqs=ZOtf(firstIndexBelow5pct+[-1 0]);
        mtfValues=mtf(firstIndexBelow5pct+[-1 0]);
        maxSpFreq=sum((1-([-1 1].*(mtfValues-mtfThreshold)./diff(mtfValues))).*spFreqs);
    else
        maxSpFreq=-ZOtf(1);
    end
    axialResolution=1/maxSpFreq;
end