%
%
%
function [deconvolvedImgs xRangeRestoration yRangeRestoration zRangeRestoration]=simulateAiryAndBesselSI()
    close all;

    recalibrate=false;
    recalculateLightSheetScan=false;

%     %     % Extra large USAF size
%     xRangeSample=single(2*[-3:.1:3]*1e-6); % propagation
%     yRangeSample=single(2*[-25:.1:25]*1e-6); % swipe
%     zRangeSample=single(2*[-25:.1:25]*1e-6); % scan

%     % USAF size
    xRangeSample=single(2*[-2.5:.1:2.5]*1e-6); % propagation
    yRangeSample=single(2*[-18.7:.1:18.7]*1e-6); % swipe
    zRangeSample=single(2*[-18.7:.1:18.7]*1e-6); % scan

%     % 'quick' calc
%     xRangeSample=single([-2.5:.1:2.5]*1e-6); % propagation
%     yRangeSample=single([-3:.1:3]*1e-6); % swipe
%     zRangeSample=single([-3:.1:3]*1e-6); % scan
    
%     alphas=[0 0   0    0  0   3 4 5 7];
%     betas= [1 .10 .05 .02 .01 1 1 1 1];
    alphas=[0 7];
    betas= [1 1];

    xOffsets=[0]*1e-6;
%     xOffsets=[0 20 25 40 50 75 80]*1e-6;
    outputDir='C:\Users\Tom\Dropbox\BesselSI\';
    outputFileName=fullfile(outputDir,'usafBeta52pctAlpha2SNR3At0and50umand100.mat');
    peakValueInRecording=0.80;

%     sampleFunctor=@getCentralPointSource; %@getUSAFX; %@getUSAFZ; % @getCentralPointSource; % @getDiskSample; %@getSphereSample; %@getTetraedronDensity; % @getSpiralDensity; % @getHomogeneousSample; @getCentralPointSource; % 
    sampleFunctor=@getUSAFX;
    
    nbParallelBeams=1;
    nbPhases=3;
    excitationWavelength=532e-9;
    fluorescenceWavelength=605e-9; % Red fluorescence spheres of Duke Scientific
    period=1e-6;
    diffractionPeriod=ceil(10e-6/period)*period;
%     besselBeta=0.02;
%     outerNA=0.42; innerNA=outerNA*(1-besselBeta); detectionNA=0.40;
    noiseToSignalRatio=.01; %10^-(15/10);
    
    config=struct();
    config.excitation=struct();
    config.excitation.objective=struct();
    config.excitation.objective.refractiveIndex=1.00;
    config.excitation.objective.numericalAperture=0.42;
%     config.excitation.objective.numericalAperture=outerNA;
    config.excitation.objective.magnification=20;
    config.excitation.objective.tubeLength=0.200;
%     config.excitation.wavelength=532e-9;
    config.excitation.wavelength=excitationWavelength;
	config.excitation.fractionOfNumericalApertureUsed=1;
    config.sample.refractiveIndex=1.40;
%     config.sample.refractiveIndex=1.33; % Fig 1C
    config.detection=struct();
    config.detection.objective=struct();
    config.detection.objective.refractiveIndex=1.00;
    config.detection.objective.numericalAperture=0.40;
%     config.detection.objective.numericalAperture=detectionNA;
    config.detection.objective.magnification=20;
    config.detection.objective.tubeLength=0.160;
    config.detection.tubeLength=0.176;
%     config.detection.wavelength=600e-9; %6.12e-7;
    config.detection.wavelength=fluorescenceWavelength;
	config.detection.fractionOfNumericalApertureUsed=1;
    config.detector=struct();
    config.detector.pixelSize=[1 1]*7.4e-6;
    
    config.detector.numberOfGrayLevels=256;
    config.detector.wellDepth=18000;
    config.sample.signalLevel=0.03;
    config.sample.backgroundLevel=0.01;
    config.detector.readOutNoise=0; %16; % electrons
    config.detector.darkCurrent=0; %200; % elec/s
    config.detector.integrationTime=30e-3; % s
        
    % Some diagnostics...
    samplePixelSize=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);    
    realMagnification=config.detection.objective.magnification*(config.detection.tubeLength/config.detection.objective.tubeLength);
    nyquistPixelPitch=realMagnification*[1 1]*0.5*config.detection.wavelength/(2*config.detection.objective.numericalAperture);
    if any(config.detector.pixelSize>nyquistPixelPitch)
        logMessage('Magnification too low! Pixel pitch (%0.1fum,%0.1fum) should be no larger than (%0.1fum,%0.1fum)!',[config.detector.pixelSize nyquistPixelPitch]*1e6);
    end
    if any(samplePixelSize(1:2)>[1 .5].*nyquistPixelPitch/realMagnification)
        logMessage('Sample space grid (%0.0fnm,%0.0fnm,%0.0fnm) too coarse, should be finer than (%0.0fnm,%0.0fnm) ',[samplePixelSize(1:3) [1 .5].*nyquistPixelPitch/realMagnification]*1e9);
    end
    
    % Simulate
    logMessage('Simulating imaging process...');
    nbLightSheetTypes=numel(alphas);
    deconvolvedImgs=cell(nbLightSheetTypes,numel(xOffsets));
    xRangeRestorations=cell(nbLightSheetTypes,numel(xOffsets));
    yRangeRestorations=cell(nbLightSheetTypes,numel(xOffsets));
    zRangeRestorations=cell(nbLightSheetTypes,numel(xOffsets));
    peakIntensities=zeros(nbLightSheetTypes,numel(xOffsets));
    meanIntensities=zeros(nbLightSheetTypes,numel(xOffsets));
    photonsPerIrradiationAtWaist=zeros(1,nbLightSheetTypes);
    for lightSheetIdx=1:nbLightSheetTypes,
        alpha=alphas(lightSheetIdx);
        besselBeta=betas(lightSheetIdx);
        for xOffsetIdx=1:numel(xOffsets),
            if xOffsetIdx>1,
                config.sample.photonsPerIrradiation=photonsPerIrradiationAtWaist(lightSheetIdx);
            else
                config.sample.photonsPerIrradiation=[]; % Will be set after the first loop
                
                if besselBeta<1
                    % Calibrate
                    cache=Cache();
                    key={config,besselBeta,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,noiseToSignalRatio};
                    if ~isfield(cache,key) || recalibrate,
                        logMessage('Calibrating Bessel SI...');
                        calibration=calibrateStructuredIllumination(config,besselBeta,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,noiseToSignalRatio,recalculateLightSheetScan);
                        cache(key)=calibration;
                    else
                        logMessage('Loading calibration from cache.');
                        calibration=cache(key);
                    end
                    componentConstantsAbsAngle=[abs(calibration.componentConstants).' angle(calibration.componentConstants).']
                end
            end
            % Simulate real sample
            if besselBeta<1
                % Simulate
                logMessage('Simulating Bessel%0.0f SR-SIM on target at %0.0fum...',[betas(lightSheetIdx)*100 xOffsets(xOffsetIdx)*1e6]);
                [deconvolvedImgs{lightSheetIdx,xOffsetIdx} xRangeRestorations{lightSheetIdx,xOffsetIdx} yRangeRestorations{lightSheetIdx,xOffsetIdx} zRangeRestorations{lightSheetIdx,xOffsetIdx} photonsPerIrradiationAtWaist(lightSheetIdx) peakIntensities(lightSheetIdx,xOffsetIdx) meanIntensities(lightSheetIdx,xOffsetIdx) irradiation(lightSheetIdx,xOffsetIdx)]=...
                    simulateStructuredIlluminationAndRestoration(config,besselBeta,calibration,sampleFunctor,xRangeSample+xOffsets(xOffsetIdx),yRangeSample,zRangeSample,xOffsets(xOffsetIdx),peakValueInRecording,recalculateLightSheetScan);
            else
                if alpha==0
                    logMessage('Simulating Conventional lightsheet on target at %0.0fum...',xOffsets(xOffsetIdx)*1e6);
                    [deconvolvedImgs{lightSheetIdx,xOffsetIdx} xRangeRestorations{lightSheetIdx,xOffsetIdx} yRangeRestorations{lightSheetIdx,xOffsetIdx} zRangeRestorations{lightSheetIdx,xOffsetIdx} photonsPerIrradiationAtWaist(lightSheetIdx) peakIntensities(lightSheetIdx,xOffsetIdx) meanIntensities(lightSheetIdx,xOffsetIdx) irradiation(lightSheetIdx,xOffsetIdx)]=...
                        simulateGaussian(config,sampleFunctor,xRangeSample+xOffsets(xOffsetIdx),yRangeSample,zRangeSample,xOffsets(xOffsetIdx),peakValueInRecording,recalculateLightSheetScan);
                else
                    logMessage('Simulating Airy%0.0f lightsheet on target at %0.0fum...',[alphas(lightSheetIdx) xOffsets(xOffsetIdx)*1e6]);
                    [deconvolvedImgs{lightSheetIdx,xOffsetIdx} xRangeRestorations{lightSheetIdx,xOffsetIdx} yRangeRestorations{lightSheetIdx,xOffsetIdx} zRangeRestorations{lightSheetIdx,xOffsetIdx} photonsPerIrradiationAtWaist(lightSheetIdx) peakIntensities(lightSheetIdx,xOffsetIdx) meanIntensities(lightSheetIdx,xOffsetIdx) irradiation(lightSheetIdx,xOffsetIdx)]=...
                        simulateAiry(config,alpha,sampleFunctor,xRangeSample+xOffsets(xOffsetIdx),yRangeSample,zRangeSample,xOffsets(xOffsetIdx),peakValueInRecording,recalculateLightSheetScan);
                end
            end
        end
        clear calibration; % Make sure not to mix up Besse SI calibrations!
    end
    
    % The following normalizes everything to first light sheet, which should be Gaussian
    totalIrradiations=[irradiation(:,:).integral]; totalIrradiations=reshape(totalIrradiations,nbLightSheetTypes,[]).'./totalIrradiations(1);
    peaksOfLightSheets=[irradiation(:).peakOfLightSheet]; peaksOfLightSheets=reshape(peaksOfLightSheets,nbLightSheetTypes,[]).'./peaksOfLightSheets(1);
    peakOfBeams=[irradiation(:).peakOfBeams]; peakOfBeams=reshape(peakOfBeams,nbLightSheetTypes,[]).'./peakOfBeams(1);
    
    [alphas; betas; totalIrradiations; peaksOfLightSheets; peakOfBeams]
    [xOffsets*1e6; peakIntensities*100; meanIntensities*100]
    
    save(outputFileName);
    
    %
    % Output
    %
    if (nargout==0)
        figs=[];
        for lightSheetIdx=1:size(deconvolvedImgs,1),
            if ~isempty(deconvolvedImgs{lightSheetIdx,1})
                figs(lightSheetIdx)=figure('Position',[50 50 1024 768],'NumberTitle','off','Name',sprintf('Lightsheet %d: alpha=%0.1f, beta=%0.0f%%.',[lightSheetIdx alphas(lightSheetIdx) betas(lightSheetIdx)*100]));
%                 [~,YOtf,ZOtf]=getSpatialFrequencyGrid(calibration.restorationPixelSize,[1 calibration.restorationSize(2:3)]);
                plotRows=1;
                plotCols=numel(xOffsets);
                axs=[]; lbs=[]; titles=[];
                for xOffsetIdx=1:numel(xOffsets),
                    deconvolvedImg=deconvolvedImgs{lightSheetIdx,xOffsetIdx};
                    xRangeRestoration=xRangeRestorations{lightSheetIdx,xOffsetIdx};
                    yRangeRestoration=yRangeRestorations{lightSheetIdx,xOffsetIdx};
                    zRangeRestoration=zRangeRestorations{lightSheetIdx,xOffsetIdx};
%                     deconvolvedImgFft=fftn(ifftshift(deconvolvedImg));
                    ySel=abs(yRangeRestoration)<=35e-6;
                    zSel=-29e-6<=zRangeRestoration & zRangeRestoration<=31e-6;
                    deconvolvedImg=deconvolvedImg(:,ySel,zSel);
                    yRangeRestoration=yRangeRestoration(ySel);
                    zRangeRestoration=zRangeRestoration(zSel);

                    xIdxAtCenter=find(abs(xRangeRestoration-xOffsets(xOffsetIdx))==min(abs(xRangeRestoration-xOffsets(xOffsetIdx))),1);

                    imgSliceYZ=permute(deconvolvedImg(xIdxAtCenter,:,:),[3 2 1]);
                    imgSliceYZ=imgSliceYZ./max(imgSliceYZ(:));
                    imwrite(imgSliceYZ(end:-1:1,:),fullfile(outputDir,sprintf('imageSliceYZ_alpha%0.0f_beta%0.0fpct_x%0.0fum.png',[alphas(lightSheetIdx) betas(lightSheetIdx)*100 xOffsets(xOffsetIdx)*1e6])));
                    
                    axs(xOffsetIdx,1)=subplot(plotRows,plotCols,xOffsetIdx+(1-1)*plotCols);
                    showImage(imgSliceYZ,[],yRangeRestoration*1e6,zRangeRestoration*1e6);
                    lbs(xOffsetIdx,1,1)=xlabel('y [\mum]'); lbs(xOffsetIdx,1,2)=ylabel('z [\mum]');
                    titles(xOffsetIdx)=title(sprintf('x = %0.0f \\mum',xOffsets(xOffsetIdx)*1e6));

%                     axs(xOffsetIdx,2)=subplot(plotRows,plotCols,xOffsetIdx+(2-1)*plotCols);
%                     showImage(deconvolvedImg(:,:,1+floor(end/2)).',-1,xRangeRestoration*1e6,yRangeRestoration*1e6);
%                     lbs(xOffsetIdx,2,1)=xlabel('x [\mum]'); lbs(xOffsetIdx,2,2)=ylabel('y [\mum]');

%                     axs(xOffsetIdx,3)=subplot(plotRows,plotCols,xOffsetIdx+(3-1)*plotCols);
%                     showImage(fftshift(permute(log10(abs(sum(deconvolvedImgFft,1))),[2 3 1])),-1,squeeze(ZOtf(1,1,:))*1e-6,YOtf(1,:,1)*1e-6);
%                     lbs(xOffsetIdx,3,1)=xlabel('\nu_z [cycles/\mum]'); lbs(xOffsetIdx,3,2)=ylabel('\nu_y [cycles/\mum]');
                end
        
                set(axs(:,1),'XTick',[-1000:10:1000],'YTick',[-1000:10:1000]);
%                 set(axs(:,3),'XTick',[-100:1:100],'YTick',[-100:1:100]);
                axis(axs(:),'equal','tight');
                set(axs(:),'YDir','normal','TickDir','out','LineWidth',2,'FontSize',14,'FontWeight','bold','FontName','Arial');
                set([lbs(:); titles(:)],'FontSize',16,'FontWeight','bold','FontName','Arial');
            end
            plot2svg(fullfile(outputDir,sprintf('usafAlpha%0.0fBeta%0.0f.svg',[alphas(lightSheetIdx),betas(lightSheetIdx)*100])),figs(lightSheetIdx));
        end
        
        clear deconvolvedImgs;
    end

end
function calibration=calibrateStructuredIllumination(config,besselBeta,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,noiseToSignalRatio,recalculateLightSheetScan)
    calibration=struct();
    calibration.nbPhases=nbPhases;
    calibration.period=period;
    calibration.nbBeams=nbParallelBeams;
    calibration.diffractionPeriod=diffractionPeriod;
    config.modulation.alpha=0;
    config.modulation.beta=besselBeta;
    peakValueInRecording=[];
    [img, xRangeDetectorInSample, yRangeDetectorInSample, zRangeDetectorInSample, photonsPerIrradiation, detectionImpulseResponses, lightSheetAtWaist, lightSheet, irradiation]=simulateGenericLightSheetScan(config,@getCentralPointSource,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,0,false,peakValueInRecording,recalculateLightSheetScan);
    clear lightSheetAtWaist lightSheet sample;
            
    clear xRangeSample yRangeSample zRangeSample;
    detectorPixelInSampleSize=voxelSizeFromRanges(xRangeDetectorInSample,yRangeDetectorInSample,zRangeDetectorInSample);
    detectorPixelInSampleSize=cast(detectorPixelInSampleSize,class(img));
    
    [~,recordingSize]=voxelSizeFromRanges(xRangeDetectorInSample,yRangeDetectorInSample,zRangeDetectorInSample);
    
    components=splitImageSetInComponents(img);
    
    % Normalize the detection impulse responses and calculate their transfer functions
    detectionImpulseResponses=ifftshift(ifftshift(ifftshift(detectionImpulseResponses,3),2),1);
    detectionImpulseResponses=detectionImpulseResponses./repmat(sum(sum(sum(detectionImpulseResponses))),[recordingSize(1:3) 1]);
    normalizedDetectionTransferFunctions=fft(fft(fft(detectionImpulseResponses,[],3),[],2),[],1);
    
    %
    % Estimate parameters
    %
    orders=[0:floor(calibration.nbPhases/2), -floor(calibration.nbPhases/2):-1];
    shiftedComponents=subPixelShift(components,yRangeDetectorInSample,orders/period);
    calibration.componentConstants=zeros([1 length(orders)]);
    calibration.componentConstants(1)=1;
    for orderIdx=2:length(orders),
        weightedComponent=shiftedComponents(:,:,:,orderIdx).*normalizedDetectionTransferFunctions(:,:,:,1);
        weightedCentralComponent=shiftedComponents(:,:,:,1).*subPixelShift(normalizedDetectionTransferFunctions(:,:,:,orderIdx),yRangeDetectorInSample,orders(orderIdx)/period);
        calibration.componentConstants(orderIdx)=sum(conj(weightedCentralComponent(:)).*weightedComponent(:))/max(eps(weightedCentralComponent(1)),sum(abs(weightedCentralComponent(:)).^2));
    end
%     componentConstantsAbs=(abs(calibration.componentConstants)+abs(calibration.componentConstants([1 end:-1:2])))/2;
%     componentConstantsArg=(angle(calibration.componentConstants)-angle(calibration.componentConstants([1 end:-1:2])))/2;
%     calibration.componentConstants=calibration.componentConstantsAbs.*exp(1i*(angle(calibration.componentConstants)-calibration.componentConstantsArg));
    
    % Increase sampling density in y and z as needed
    calibration.restorationUpSamplingFactor=[1 2 1];
    calibration.restorationPixelSize=detectorPixelInSampleSize./calibration.restorationUpSamplingFactor;
    calibration.restorationSize=recordingSize.*calibration.restorationUpSamplingFactor;
    normalizedDetectionTransferFunctions=circshift(normalizedDetectionTransferFunctions,[floor(recordingSize./2) 0]);
    normalizedDetectionTransferFunctions(calibration.restorationSize(1),calibration.restorationSize(2),calibration.restorationSize(3),1)=0;
    normalizedDetectionTransferFunctions=circshift(normalizedDetectionTransferFunctions,-[floor(recordingSize./2) 0]);
    detectionImpulseResponses=ifft(ifft(ifft(normalizedDetectionTransferFunctions,[],3),[],2),[],1);
    yRangeRestoration=yRangeDetectorInSample(1+floor(end/2))+([1:calibration.restorationSize(2)]-1-floor(calibration.restorationSize(2)/2))*calibration.restorationPixelSize(2);
    
    %
    % Construct filter terms for each order
    %  
    lateralCutOffSpatialFrequency=(2*config.detection.objective.numericalAperture)/config.detection.wavelength;
    axialCutOffSpatialFrequencyExc=(2*config.excitation.objective.numericalAperture)/config.excitation.wavelength;
    axialCutOffSpatialFrequencyDet=(2*config.sample.refractiveIndex*(1-cos(asin(config.detection.objective.numericalAperture/config.sample.refractiveIndex))))/config.detection.wavelength;
    cutOffSpatialFrequencies=[lateralCutOffSpatialFrequency lateralCutOffSpatialFrequency+max(abs(orders))/period axialCutOffSpatialFrequencyExc+axialCutOffSpatialFrequencyDet];
%     cutOffSpatialFrequencies=[lateralCutOffSpatialFrequency [1 1]*lateralCutOffSpatialFrequency+max(abs(orders))/period];
    [XOtf,YOtf,ZOtf]=getSpatialFrequencyGrid(calibration.restorationPixelSize,calibration.restorationSize(1:3));
    calibration.filterTerms=zeros([calibration.restorationSize(1:3) length(orders)],class(detectionImpulseResponses));
    for orderIdx=1:length(orders),
        detectionTransferFunctions=zeros(size(calibration.filterTerms),class(detectionImpulseResponses));
        for orderPrimeIdx=1:length(orders),
            deflection=exp(2i*pi*(orders(orderPrimeIdx)-orders(orderIdx))*yRangeRestoration/period);
            deflection=repmat(ifftshift(deflection),[calibration.restorationSize(1) 1 calibration.restorationSize(3) 1]);
            detectionTransferFunctions(:,:,:,orderPrimeIdx)=calibration.componentConstants(orderPrimeIdx)*fftn(detectionImpulseResponses(:,:,:,orderIdx).*deflection);
        end
%         clear deflection detectionImpulseResponses;
        denominator=sum(abs(detectionTransferFunctions).^2,4)+noiseToSignalRatio.^2;
%         clear detectionTransferFunctions;
        filter=conj(calibration.componentConstants(orderIdx)*normalizedDetectionTransferFunctions(:,:,:,orderIdx))./denominator;
        % Remove high frequency noise now
        fRel=sqrt(...
            (XOtf/cutOffSpatialFrequencies(1)).^2 + ...
            ((YOtf+orders(orderIdx)/period)/cutOffSpatialFrequencies(2)).^2 +...
            (ZOtf/cutOffSpatialFrequencies(3)).^2 );
        apodization=ifftshift(max(0,1-fRel));
        clear fRel;
        calibration.filterTerms(:,:,:,orderIdx)=filter.*apodization;
        clear filter apodization;
    end
    clear XOtf YOtf ZOtf;
    
end
function [deconvolvedImg xRangeRestoration yRangeRestoration zRangeRestoration photonsPerIrradiation peakIntensity meanIntensity irradiation]=simulateStructuredIlluminationAndRestoration(config,besselBeta,calibration,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,xOffset,peakValueInRecording,recalculateLightSheetScan)
    config.modulation.alpha=0;
    config.modulation.beta=besselBeta;
    [img, xRangeDetectorInSample, yRangeDetectorInSample, zRangeDetectorInSample, photonsPerIrradiation, detectionImpulseResponses, lightSheetAtWaist, lightSheet, irradiation]=simulateGenericLightSheetScan(config,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,calibration.nbPhases,calibration.period,calibration.nbBeams,calibration.diffractionPeriod,xOffset,true,peakValueInRecording,recalculateLightSheetScan);
    clear xRangeSample yRangeSample zRangeSample;
    [~,recordingSize]=voxelSizeFromRanges(xRangeDetectorInSample,yRangeDetectorInSample,zRangeDetectorInSample);
    
    peakIntensity=max(img(:));
    meanIntensity=mean(img(:));
    
    logMessage('Peak intensity of recorded image %0.3f%%.',peakIntensity*100);
    
    %
    % Reconstruct the sample image
    %    
  
    % Split image in components in frequency space, leaving 'em in place
    components=splitImageSetInComponents(img);
    
    components=circshift(components,[floor(recordingSize./2) 0]);
    components(calibration.restorationSize(1),calibration.restorationSize(2),calibration.restorationSize(3),1)=0;
    components=circshift(components,-[floor(recordingSize./2) 0]);
    
    xRangeRestoration=xRangeDetectorInSample(1+floor(end/2))+([1:calibration.restorationSize(1)]-1-floor(calibration.restorationSize(1)/2))*calibration.restorationPixelSize(1);
    yRangeRestoration=yRangeDetectorInSample(1+floor(end/2))+([1:calibration.restorationSize(2)]-1-floor(calibration.restorationSize(2)/2))*calibration.restorationPixelSize(2);
    zRangeRestoration=zRangeDetectorInSample(1+floor(end/2))+([1:calibration.restorationSize(3)]-1-floor(calibration.restorationSize(3)/2))*calibration.restorationPixelSize(3);
    clear detectorPixelInSampleSize recordingSize xRangeDetectorInSample yRangeDetectorInSample zRangeDetectorInSample;
    
    %
    % Generalized Wiener filter application
    %
    % Restore components one by one in situ
    filteredComponents=components.*calibration.filterTerms;
    % Shift components back to original spot
    orders=[0:floor(calibration.nbPhases/2), -floor(calibration.nbPhases/2):-1];
    filteredComponents=subPixelShift(filteredComponents,yRangeRestoration,orders/calibration.period);
    % Combine
    deconvolvedImgFft=sum(filteredComponents,4);
    deconvolvedImg=fftshift(ifftn(deconvolvedImgFft,'symmetric'));
end

function [img xRangeDetectorInSample yRangeDetectorInSample zRangeDetectorInSample photonsPerIrradiation peakIntensity meanIntensity irradiation]=simulateGaussian(config,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,xOffset,peakValueInRecording,recalculateLightSheetScan)
    % Simulate recording
    config.modulation.alpha=0;
    config.modulation.beta=1;
    simulateNoise=true;
    [img, xRangeDetectorInSample, yRangeDetectorInSample, zRangeDetectorInSample, photonsPerIrradiation, detectionImpulseResponses, lightSheetAtWaist, lightSheet, irradiation]=...
        simulateGenericLightSheetScan(config,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,1,0,1,0,xOffset,simulateNoise,peakValueInRecording,recalculateLightSheetScan);
     
    peakIntensity=max(img(:));
    meanIntensity=mean(img(:));
   
end

function [deconvolvedImg xRangeRestoration yRangeRestoration zRangeRestoration photonsPerIrradiation peakIntensity meanIntensity irradiation]=simulateAiry(config,alpha,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,xOffset,peakValueInRecording,recalculateLightSheetScan)
    % Simulate recording
    config.modulation.alpha=alpha;
    config.modulation.beta=1;
    simulateNoise=true;
    [img, xRangeDetectorInSample, yRangeDetectorInSample, zRangeDetectorInSample, photonsPerIrradiation, detectionImpulseResponses, lightSheetAtWaist, lightSheet, irradiation]=...
        simulateGenericLightSheetScan(config,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,1,0,1,0,xOffset,simulateNoise,peakValueInRecording,recalculateLightSheetScan);
     
    peakIntensity=max(img(:));
    meanIntensity=mean(img(:));
    
    
%     realMagnification=config.detection.objective.magnification*(config.detection.tubeLength/config.detection.objective.tubeLength);
%     lightSheetDownSampled=resampleImage(lightSheet,xRangeSample,yRangeSample,config.detector.pixelSize/realMagnification);
%     lightSheetDownSampled=sum(lightSheetDownSampled,2);
%     lightSheetDownSampled=lightSheetDownSampled./max(lightSheetDownSampled(:));
    
    % Simulate deconvolution
    config.stagePositions.target=zRangeSample;
    config.detector.center=-[xRangeDetectorInSample(1+floor(end/2)) yRangeDetectorInSample(1+floor(end/2))];
    
    % Swap (x,y,z)=(propagation,swipe,scan) as in paper to (swipe,propagation,scan) as in processing code
    img=permute(img,[2 1 3]);
    config.detector.center(1:2)=config.detector.center([2 1]);
    [deconvolvedImg lightSheetDeconvFilter lightSheetOtf ZOtf yRangeRestoration,xRangeRestoration,zRangeRestoration tRange lightSheetPsf]=deconvolveRecordedImageStack(img,config);
    deconvolvedImg=ipermute(deconvolvedImg,[2 1 3]);
    
%     lightSheetPsf=ipermute(lightSheetPsf,[2 1 3]);
end


%
% shiftedData=subPixelShift(data,yRange,frequencies)
%
% Shifts a matrix via its Fourier transform
function shiftedData=subPixelShift(data,yRange,frequencies)
    dataSize=[size(data,1), size(data,2), size(data,3)];
    nbOrders=numel(frequencies);
    deflection=repmat(ifftshift(yRange),[1 1 1 nbOrders]).*repmat(permute(frequencies,[1 4 3 2]),[1 length(yRange) 1 1]);
    deflection=repmat(exp(2i*pi*deflection),[dataSize(1) 1 dataSize(3) 1]);
    dataIFft=ifft(ifft(ifft(data,[],3),[],2),[],1);
    shiftedData=fft(fft(fft(dataIFft.*deflection,[],3),[],2),[],1);
end
    
function sample=getSpiralDensity(xRangeSample,yRangeSample,zRangeSample)
    diameter=20e-6;
    period=2*diameter;
    tubeDiameter=10e-6;
    thickness=1e-6;
    
    nbSpirals=3;
    
    innerR2=abs(tubeDiameter/2-thickness/2)^2;
    outerR2=abs(tubeDiameter/2+thickness/2)^2;
    
    % Make all horizontal
    xRangeSample=xRangeSample(:).'; yRangeSample=yRangeSample(:).'; zRangeSample=zRangeSample(:).';
    
    [~,outputSize]=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);
    
    sample=false(outputSize);
    
    for xIdx=1:outputSize(1),
        thetas=2*pi*(xRangeSample(xIdx)/period+([1:nbSpirals].'-1)./nbSpirals);
        spiralPos=[thetas*0 cos(thetas) sin(thetas)]*diameter/2;
        inAnySpiral=false(outputSize(2:3));
        for spiralIdx=1:size(spiralPos,1),
            R2=repmat(abs(yRangeSample-spiralPos(spiralIdx,2)).^2.',[1 outputSize(3)])+repmat(abs(zRangeSample-spiralPos(spiralIdx,3)).^2,[outputSize(2) 1]);
            inSpiral=R2>=innerR2 & R2<=outerR2;
            inAnySpiral=inAnySpiral|inSpiral;
        end
        sample(xIdx,:,:)=inAnySpiral;
    end
    
end

function sample=getUSAFX(xRangeSample,yRangeSample,zRangeSample)
    sample=getUSAFZ(zRangeSample,yRangeSample,xRangeSample);
    sample=permute(sample,[3 2 1]);
end

function sample=getUSAFZ(xRangeSample,yRangeSample,zRangeSample)
    sample=zeros([numel(xRangeSample) numel(yRangeSample) numel(zRangeSample)]);
    img=1-getTestImage('usaf1951_375x375.png').';
    img=img(:,end:-1:1).';
    % 75mm/375px => 200um/px, largest line (-2,1) 10x2mm
    % 37.5um/375px => 100nm/px, largest line (-2,1) 5x1um
    % Crop or zero pad image to fit the sample size
    sampleSize=size(sample);
%     % Down scale the test image first 
%     while size(img,1)>sampleSize(1),
%         img=0.5*(img(1:2:end-1,:)+img(2:2:end,:));
%     end
%     while size(img,2)>sampleSize(2),
%         img=0.5*(img(:,1:2:end-1)+img(:,2:2:end));
%     end
    topLeftDiff=floor(sampleSize(1:2)/2)-floor(size(img)/2);
    if size(img,1)<size(sample,1)
        img(size(sample,1),1)=0;
        img=circshift(img,[topLeftDiff(1) 0]);
    else
        img=circshift(img,[topLeftDiff(1) 0]);
        img=img(1:size(sample,1),:);
    end
    if size(img,2)<size(sample,2)
        img(1,size(sample,2))=0;
        img=circshift(img,[0 topLeftDiff(2)]);
    else
        img=circshift(img,[0 topLeftDiff(2)]);
        img=img(:,1:size(sample,2));
    end
    
    sample=repmat(img,[1 1 numel(zRangeSample)]);
end

function sample=getHomogeneousSample(xRangeSample,yRangeSample,zRangeSample)
    sample=ones([numel(xRangeSample) numel(yRangeSample) numel(zRangeSample)]);
end

function sample=getCentralPointSource(xRangeSample,yRangeSample,zRangeSample)
    sample=zeros([numel(xRangeSample) numel(yRangeSample) numel(zRangeSample)]);
    sample(abs(xRangeSample)==min(abs(xRangeSample)),abs(yRangeSample)==min(abs(yRangeSample)),abs(zRangeSample)==min(abs(zRangeSample)))=1;
    
%     sample(abs(xRangeSample-1.5e-6)==min(abs(xRangeSample-1.5e-6)),abs(yRangeSample)==min(abs(yRangeSample)),abs(zRangeSample)==min(abs(zRangeSample)))=1;
%     sample(abs(xRangeSample)==min(abs(xRangeSample)),abs(yRangeSample-1e-6)==min(abs(yRangeSample-1e-6)),abs(zRangeSample)==min(abs(zRangeSample)))=1;
%     sample(abs(xRangeSample)==min(abs(xRangeSample)),abs(yRangeSample)==min(abs(yRangeSample)),abs(zRangeSample-2e-6)==min(abs(zRangeSample-2e-6)))=1;
end

function sample=getDiskSample(xRangeSample,yRangeSample,zRangeSample)
    sample=getSphereSample(xRangeSample,yRangeSample,zRangeSample);
    sample(:,:,[1:floor(end/2), 2+floor(end/2):end])=0; % Blacken out everything but z==0
end

function sample=getSphereSample(xRangeSample,yRangeSample,zRangeSample)
    diameter=20e-6;
    thickness=3e-6;

    % Ensure all horizontal
    xRangeSample=xRangeSample(:).'; yRangeSample=yRangeSample(:).'; zRangeSample=zRangeSample(:).';
    
    [~, outputSize]=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);
    
    innerR2=abs(diameter/2-thickness/2)^2;
    outerR2=abs(diameter/2+thickness/2)^2;
    
    sample=false(outputSize);
    
    for xIdx=1:outputSize(1),
        R2=abs(xRangeSample(xIdx)).^2+...
            repmat(abs(yRangeSample).^2.',[1 outputSize(3)])+repmat(abs(zRangeSample).^2,[outputSize(2) 1]);
        inShell=R2>=innerR2 & R2<=outerR2;
        sample(xIdx,:,:)=inShell;
    end
end

function sample=getTetraedronDensity(xRangeSample,yRangeSample,zRangeSample)
    diameter=20e-6;
    thickness=1e-6;
    centerOffset=[0 -5e-6 -2.5e-6];
    
    % Ensure all horizontal
    xRangeSample=xRangeSample(:).'; yRangeSample=yRangeSample(:).'; zRangeSample=zRangeSample(:).';
    
    [~, outputSize]=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);
    
    innerR2=abs(diameter/2-thickness/2)^2;
    outerR2=abs(diameter/2+thickness/2)^2;
    
    sample=false(outputSize);
    
    spherePos=[0 -1/sqrt(6) 2/sqrt(3); -1 -1/sqrt(6) -1/sqrt(3); 1 -1/sqrt(6) -1/sqrt(3); 0 sqrt(3/2) 0]*diameter/2;
    spherePos=spherePos+repmat(centerOffset,[size(spherePos,1) 1]);
    
    for xIdx=1:outputSize(1),
        inAnyShell=false(outputSize(2:3));
        for sphereIdx=1:size(spherePos,1),
            R2=abs(xRangeSample(xIdx)-spherePos(sphereIdx,1)).^2+...
                repmat(abs(yRangeSample-spherePos(sphereIdx,2)).^2.',[1 outputSize(3)])+repmat(abs(zRangeSample-spherePos(sphereIdx,3)).^2,[outputSize(2) 1]);
            inShell=R2>=innerR2 & R2<=outerR2;
            inAnyShell=inAnyShell|inShell;
        end
        sample(xIdx,:,:)=inAnyShell;
    end
    
end

% the object is centered at xOffset
function [img, xRangeDetectorInSample, yRangeDetectorInSample, zRangeDetectorInSample, photonsPerIrradiation, detectionImpulseResponses, lightSheetAtWaist, lightSheet, irradiation]=simulateGenericLightSheetScan(config,sampleFunctor,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,xOffset,simulateNoise,peakValueInRecording,recalculate)
    % Load synthetic sample
    sample=sampleFunctor(xRangeSample-xOffset,yRangeSample,zRangeSample);
    sample=cast(sample,class(xRangeSample));
    [~, sampleSize]=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);
    
    % Check if we can get this from the cache
    key={config,sample,xRangeSample,yRangeSample,zRangeSample,nbPhases,period,nbParallelBeams,diffractionPeriod,xOffset};
    cache=Cache();
    if ~isfield(cache,key) || recalculate
        %
        % Calculate light sheet
        %
        pupilFunctor=@(U,V) (U.^2+V.^2>=(1-config.modulation.beta)^2).*exp(2i*pi*config.modulation.alpha*(U.^3+V.^3));
        lightSheet=calcLightSheet(config,xRangeSample,yRangeSample,zRangeSample,pupilFunctor,diffractionPeriod,nbParallelBeams,period,nbPhases);
        lightSheetAtWaist=calcLightSheet(config,[0 max(abs(xRangeSample))],yRangeSample,zRangeSample,pupilFunctor,diffractionPeriod,nbParallelBeams,period,1); % Force the same pupil grid size by using the same maximum defocus
        lightSheetAtWaist=lightSheetAtWaist(1,:,:); % Calculate the light sheet at the waist position in case xRangeSample does not contain 0
        lightSheet=lightSheet./max(lightSheetAtWaist(:));  % Maximum excitation at waist
        lightSheetAtWaist=lightSheetAtWaist./max(lightSheetAtWaist(:));
    %     % Use plain regular structured illumination instead:
    %     lightSheet=zeros([sampleSize nbPhases],class(xRangeSample));
    %     for idx=1:nbPhases
    %         lightSheet(:,:,:,idx)=repmat(0.5+0.5*cos(2*pi*((idx-1)/nbPhases+yRangeSample/period)),[sampleSize(1) 1 sampleSize(3)]);
    %     end

        %
        % Calculate detection PSF
        %
        logMessage('Calculating detection PSF...');
        centerCoord=[xRangeSample(1+floor(end/2)) yRangeSample(1+floor(end/2)) zRangeSample(1+floor(end/2))];
        detectionPsf=calcVectorialPsf(xRangeSample-centerCoord(1),yRangeSample-centerCoord(2),zRangeSample-centerCoord(3),config.detection.wavelength,1/sqrt(2),1i/sqrt(2),...
                        config.detection.objective.numericalAperture,config.sample.refractiveIndex,...
                        config.detection.objective.magnification,config.detection.objective.tubeLength);
        detectionOtf2D=fft2(ifftshift(ifftshift(detectionPsf,1),2)); % Unity energy propagating

        %
        % Create the set of light sheet fluorescence projections
        %
        logMessage('Starting light sheet recording in z-interval [%0.3fum %0.3fum]...',[zRangeSample([1 end])*1e6]);
        stopWatch=tic();
        img=zeros([sampleSize(1:3) nbPhases],class(sample));
        for zIdx=1:sampleSize(3),
            for phaseIdx=1:nbPhases,
                fluorescence=sample(:,:,max(1,min(end,[1:end]+(zIdx-1)-floor(end/2)))).*lightSheet(:,:,:,phaseIdx);
                imgSlice=sum(ifft2(fft2(fluorescence).*detectionOtf2D,'symmetric'),3);
                img(:,:,zIdx,phaseIdx)=imgSlice;
            end
            if any(zIdx==[1 floor(sampleSize(3)/3) floor(sampleSize(3)*2/3)])
                secondsLeft=(toc(stopWatch)/zIdx)*(sampleSize(3)-zIdx);
                timeLeft=sprintf('%2.0f:%02.0f:%02.0f',[floor(secondsLeft/60^2) floor(mod(secondsLeft,60^2)/60) mod(secondsLeft,60)]);
                finishTime=addtodate(now(),round(secondsLeft*1000),'millisecond');
                logMessage(['Estimated finishing time: ',datestr(finishTime,'HH:MM:SS'),', time left:',timeLeft,'.']);
            end
        end
        clear detectionOtf2D;
        cache(key)={img,detectionPsf,lightSheet,lightSheetAtWaist};
    else
        logMessage('Loading image from cache...');
        values=cache(key);
        [img,detectionPsf,lightSheet,lightSheetAtWaist]=values{:};
        clear values;
    end
    
    %
    % Simulate the detection
    %
    logMessage('Simulating the detection...');
    realMagnification=config.detection.objective.magnification*(config.detection.tubeLength/config.detection.objective.tubeLength);
    [img xRangeDetectorInSample yRangeDetectorInSample]=resampleImage(img,xRangeSample,yRangeSample,config.detector.pixelSize/realMagnification);
    %Scale the image intensity as requested
    if isfield(config,'sample') && isfield(config.sample,'photonsPerIrradiation') && ~isempty(config.sample.photonsPerIrradiation)
        photonsPerIrradiation=config.sample.photonsPerIrradiation;
        peakValueInRecording=[];
    else
        if ~isempty(peakValueInRecording)
            photonsPerIrradiation=peakValueInRecording*config.detector.wellDepth/max(eps(class(img)),max(img(:)));
        else
            photonsPerIrradiation=1;
        end
    end
    img=img*photonsPerIrradiation;
    lightSheetAtWaist=lightSheetAtWaist*photonsPerIrradiation;
    lightSheet=lightSheet*photonsPerIrradiation;
    
    % Store some irradiation characteristics
    irradiation=struct();
    irradiation.integral=sum(lightSheet(:));
    irradiation.peakOfLightSheet=max(lightSheet(:));
    if period>0
        nbBeamsPerSwipe=(diff(yRangeSample([1 end]))+diff(yRangeSample(1:2)))/period;
        irradiation.peakOfBeams=irradiation.peakOfLightSheet*nbBeamsPerSwipe/nbParallelBeams;
    else
        irradiation.peakOfBeams=irradiation.peakOfLightSheet;
    end
    
    zRangeDetectorInSample=zRangeSample;
    if simulateNoise
        % Add noise
        darkNoiseElectrons=config.detector.readOutNoise+config.detector.darkCurrent*config.detector.integrationTime;
        outputClass=class(img);
        img=poissrnd(double(max(0,img+darkNoiseElectrons))); % Simulate photon noise, under sampling may cause negative values
        img=cast(img,outputClass);
        if any(img(:)>config.detector.wellDepth)
            logMessage('%.0f pixels saturated! Peak value %f/%f=%0.3f',[sum(img(:)>config.detector.wellDepth) max(img(:)) config.detector.wellDepth max(img(:))/config.detector.wellDepth]);
            img=min(config.detector.wellDepth,img); % Simulate detector saturation
        end
        img=floor(img*(config.detector.numberOfGrayLevels-1)/config.detector.wellDepth)/(config.detector.numberOfGrayLevels-1); % Simulate ADC
    end
    
    %
    % Calculate the detection transfer functions
    %
    logMessage('Calculating the transfer functions...');
    orders=[0:floor(nbPhases/2), -floor(nbPhases/2):-1];
    lightSheetAtWaistFft=fftn(ifftshift(lightSheetAtWaist));
    axialIlluminationsFft=subPixelShift(repmat(lightSheetAtWaistFft,[1 1 1 numel(orders)]),yRangeSample,orders/period);
    axialIlluminationsFft=axialIlluminationsFft(1,1,:,:);
    axialIlluminations=fftshift(ifft(axialIlluminationsFft(1,1,:,:),[],3,'symmetric'),3);
    detectionImpulseResponses=repmat(detectionPsf,[1 1 1 numel(orders)]).*repmat(axialIlluminations,[sampleSize(1:2) 1 1]);
    detectionImpulseResponses=resampleImage(detectionImpulseResponses,xRangeSample,yRangeSample,config.detector.pixelSize/realMagnification);
end

% Simulate detection by a pixel array, taking into account the square shape of pixels
function [resampledImg xRangeDetectorInSample yRangeDetectorInSample pixelTransferFunction]=resampleImage(img,xRangeSample,yRangeSample,detectorPixelPitchInSample)
    imgFft2D=fft2(img); clear img;
    [resampledImg xRangeDetectorInSample yRangeDetectorInSample pixelTransferFunction]=resampleImageFromFft2D(imgFft2D,xRangeSample,yRangeSample,detectorPixelPitchInSample);
end
function [resampledImg xRangeDetectorInSample yRangeDetectorInSample pixelTransferFunction]=resampleImageFromFft2D(imgFft2D,xRangeSample,yRangeSample,detectorPixelPitchInSample)
    dataSize=size(imgFft2D);
    if numel(dataSize)<4
        dataSize(4)=1;
    end
    nbSlices=prod(dataSize(3:end));
    
    % Resample
    samplePitch=[diff(xRangeSample(1:2)) diff(yRangeSample(1:2))];
    detectorGridSize=ceil(dataSize(1:2).*samplePitch./detectorPixelPitchInSample);
    [XOtfSample,YOtfSample]=getSpatialFrequencyGrid(samplePitch,dataSize(1:2));
    [XOtfDetector,YOtfDetector]=getSpatialFrequencyGrid(detectorPixelPitchInSample,detectorGridSize);
    
    % Calculate the transfer function for a square pixel of 100% fill factor
    pixelTransferFunction=sinc(XOtfDetector.*detectorPixelPitchInSample(1)).*sinc(YOtfDetector.*detectorPixelPitchInSample(2));
    pixelTransferFunction=0*pixelTransferFunction+1; % Assume 'perfect' pixels
    
    % Assume no aliassing (no replicas)
    resampledImg=zeros([detectorGridSize dataSize(3:end)],class(imgFft2D));
    for sliceIdx=1:nbSlices,
        imgFft=fftshift(imgFft2D(:,:,sliceIdx));
        resampledImgFft=interpn(XOtfSample,YOtfSample,imgFft,XOtfDetector,YOtfDetector,'*cubic',0);
        pixelatedImgFft=resampledImgFft.*pixelTransferFunction;
        resampledImg(:,:,sliceIdx)=ifft2(ifftshift(pixelatedImgFft),'symmetric');
    end
    
    if nargout>1
        xRangeDetectorInSample=xRangeSample(1+floor(end/2))+([1:detectorGridSize(1)]-1-floor(detectorGridSize(1)/2)).*detectorPixelPitchInSample(1);
        yRangeDetectorInSample=yRangeSample(1+floor(end/2))+([1:detectorGridSize(2)]-1-floor(detectorGridSize(2)/2)).*detectorPixelPitchInSample(2);
    end
end

% If period==0, continuous light sheet illumination is used
function lightSheet=calcLightSheet(config,xRangeSample,yRangeSample,zRangeSample,pupilFunctor,diffractionPeriod,nbParallelBeams,period,nbPhases)
    if (nbParallelBeams<=1)
        nbParallelBeams=1;
        diffractionPeriod=2*diff(yRangeSample([1 end])); % Fill the whole data cube with beams
    end
    structuredIllumination=period>0;
    [~, dataSize]=voxelSizeFromRanges(xRangeSample,yRangeSample,zRangeSample);

    [~,psfField]=calcVectorialPsf(yRangeSample,zRangeSample,xRangeSample,config.excitation.wavelength,@(U,V) 1/sqrt(2)*pupilFunctor(U,V),@(U,V) 1i/sqrt(2)*pupilFunctor(U,V),...
                    config.excitation.objective.numericalAperture,config.sample.refractiveIndex,...
                    config.excitation.objective.magnification,config.excitation.objective.tubeLength);
	psfField=ipermute(psfField,[2 3 1 4]);
    psfField(1,end*2,1,1)=0; % Zero pad

    if nbParallelBeams>1
        psfFieldFftY=fft(psfField,[],2);
        clear psfField;

        % Create diffractive deflection of nbParallelBeams beams
        logMessage('Passing beam through diffraction grating...');
        gratingFft=0*yRangeSample;
        gratingFft(end*2)=0;
        for diffBeamPos=diffractionPeriod*([1:nbParallelBeams]-floor(1+nbParallelBeams/2)),
            if abs(diffBeamPos)<-2*yRangeSample(1),
                gratingFft=gratingFft+exp(2i*pi*(([1:2*length(yRangeSample)]-1-length(yRangeSample))/(2*length(yRangeSample)))*diffBeamPos/diff(yRangeSample(1:2)));
            end
        end
        gratingFft=ifftshift(gratingFft)./sqrt(nbParallelBeams);

        psf=ifft(psfFieldFftY.*repmat(gratingFft(:).',[dataSize(1) 1 dataSize(3) 3]),[],2);
        clear psfFieldFftY gratingFft;
    else
        psf=psfField;
		clear psfField;
    end
    psf=sum(abs(psf).^2,4);
    
    if structuredIllumination
        % Create intensity structured illumination and phase shifts
        logMessage('Patterning light sheet intensity...');
    else
        logMessage('Swiping beam into a light sheet.');
    end
    psfFft=fft(psf,[],2);
    clear psf;
    lightSheetFft=zeros([dataSize.*[1 2 1] nbPhases],class(psfFft));
    for phaseIdx=1:nbPhases, % Generate phases
        phaseOffset=(phaseIdx-1)/nbPhases;
        
        combFft=zeros(1,2*dataSize(2));
        if structuredIllumination
            beamPositions=period*(phaseOffset+[1:round(diffractionPeriod/period)]-1-floor(round(diffractionPeriod/period)/2));
            for beamPos=beamPositions,
                combFft=combFft+exp(2i*pi*(([1:2*dataSize(2)]-1-dataSize(2))/(2*dataSize(2)))*beamPos/diff(yRangeSample(1:2)));
            end
            combFft=ifftshift(combFft)./length(beamPositions);
        else
            combFft(1)=dataSize(2);
        end
        
        lightSheetFft(:,:,:,phaseIdx)=psfFft.*repmat(combFft(:).',[dataSize(1) 1 dataSize(3)]);
    end
    clear psfFft;
    lightSheet=ifft(lightSheetFft,[],2,'symmetric');
    clear lightSheetFft;
            
    lightSheet=lightSheet(:,1:end/2,:,:); % Crop again
end

%
% [voxelSize, gridSize]=voxelSizeFromRanges(varargin)
% 
% Input arguments: a set of vectors, one for each dimension of uniformely increasing coordinates
% Output: two vectors containing:
%     voxelSize: the dimensions of each voxel
%      gridSize: the number of voxels per dimension
%
%
function [voxelSize, gridSize]=voxelSizeFromRanges(varargin)
    if nargin>0
        voxelSize=zeros(1,numel(varargin),class(varargin{1}));
        gridSize=zeros(1,numel(varargin));
        for dimIdx=1:numel(varargin)
            rng=varargin{dimIdx};
            switch numel(rng)
                case 0
                    % Empty, ignore everything after this element
                    voxelSize=voxelSize(1:dimIdx-1);
                    gridSize=gridSize(1:dimIdx-1);
                    return;
                case 1
                    voxelSize(dimIdx)=1;
                otherwise
                    voxelSize(dimIdx)=diff(rng(1:2));
            end
            gridSize(dimIdx)=numel(rng);
        end
    else
        voxelSize=[];
        gridSize=[];
    end
end

% Splits a phase shifted illumination image set into its various Fourier
% components without recentering the components
function components=splitImageSetInComponents(img)
    recordingSize=[size(img,1) size(img,2) size(img,3)];
    nbPhases=size(img,4);
    orders=[0:floor(nbPhases/2), -floor(nbPhases/2):-1];
    components=zeros([recordingSize(1:3) length(orders)],class(img));
    for orderIdx=1:length(orders),
        phaseModulation=exp(2i*pi*(([1:nbPhases]-1)/nbPhases)*orders(orderIdx));
        phaseModulation=cast(phaseModulation,class(img));
        phaseModulation=repmat(permute(phaseModulation,[1 4 3 2]),[recordingSize(1:3) 1]);
        % Modulate to select the order of interest
        component=fft(fft(fft(ifftshift(img.*phaseModulation),[],3),[],2),[],1); % Don't transform fourth dimension
        components(:,:,:,orderIdx)=sum(component,4); % Integrate to remove interference from other orders
        clear component;
    end
end