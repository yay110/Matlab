%
% fitLightSheetModelToRecording('C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-10_focus_on_left cylindrial lens at 48.5mm\Projections of 2014-06-10_11-29-05_reflected stepSize185nm.png');
% fitLightSheetModelToRecording('C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-10_focus_on_left cylindrial lens at 49mm\Projections of 2014-06-10_13-00-01_reflected stepSize185nm.png');
%
function fitLightSheetModelToRecording(lightSheetFileName,stepSize)
    close all;
    if nargin<1,
        %lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\Projections of 2014-05-29_17-31-55_reflected__stepsize500nm_pixelPitch185nm.png';
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\mean projection stepSize500nm.png';
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-10_focus_on_left cylindrial lens at 48.5mm\Projections of 2014-06-10_11-29-05_reflected stepSize185nm.png';
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-10_focus_on_left cylindrial lens at 49mm\Projections of 2014-06-10_13-00-01_reflected stepSize185nm.png';
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\10X objective 47.5mm\Projections of 2014-06-11_16-24-51_reflected stepSize740nm.png'
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-12 40X\light sheet stepSize185nm.png';
%         lightSheetFileName='C:\Users\Tom\Dropbox\OpenSPIM\data\2014-6-13 New layout 40X\2014-06-13_12-43-41_beam stepSize185nm -2.png';
%         lightSheetFileName='C:\Users\Zhengyi\Dropbox\OpenSPIM\data\2014-06-17 even illumination 40X\Projections of 2014-06-17_11-44-19_even stepSize185nm.png';
%         lightSheetFileName='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\2014-07-03\2014-07-03_13-04-11_8 beam stepSize185nm-croped.png';
%         lightSheetFileName='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\2014-06-19_Cell data\2014-06-19_20-14-50_beam stepSize185nm.png';
%lightSheetFileName='G:\2015-03-20 Airy beads for de-con\Projections of 13-00-50 stepSize325nm.png';
lightSheetFileName='F:\2015-12-07 Decon\2.6mm\AVG_2015-12-10_11-07-35 mirror stepSize325nm.png';
    end
    if nargin<2 || isempty(stepSize),
        tkns=regexpi(lightSheetFileName,'stepSize\s?([\d\.]+)\s?nm','tokens');
        stepSize=str2double(tkns{1}{1})*1e-9;
        %stepSize=stepSize*1.03; % TODO remove
        clear tkns;
    end
    
    outputFileName=strcat(lightSheetFileName(1:end-4),'_optimizationValues5_phaseAndAmplitude.mat');
    
    %To define which optimization to run.
    %For new data, best to run global first to find a starting point. 
    globalOptimization = true;
    localOptimization = false;
    globalOptimizationMinutes = 2;
    %reset starting potint (true) 
    %or use the best point from last optimization (false)
    reprocess=false;
    
    nbParamsToUse=4;
    
    %choose three Z position for fitting
    zSectionPositions=[-50:40:50]*1e-6;
    %take data from avarage of total length 2.5um.
    sectionWidth=2.5e-6;
    
    if exist(outputFileName) && ~reprocess,
        load(outputFileName,'params','err');
        params0=params;
        clear params;
    else
        params0=[0,0.8,0,5,0,0,0,0]; %  background, openFraction, defocus, alpha, fourth, fifth, amplitudeTilt, ...
    end
    
    if numel(params0)<9
        params0(9)=0;
    end
    parameterDescriptions = {'background', [0 .20], 0.001, params0(1)
                    'openFraction', [0.7 1.00], 0.01, params0(2)
                    'defocus', [-50 50], 0.1, params0(3)
                    'alpha', [5 10], 0.01, params0(4)
                    'fourth', [-10 10], 0.01, params0(5)
                    'fifth', [-10 10], 0.01, params0(6)
                    'amplitudeTilt', [-1 1], 0.01, params0(7)
                    'amplitudeQuadratic', [-1 1], 0.01, params0(8)
                    'amplitudeCubic', [-1 1], 0.01, params0(9)
                    };
	params0=params0(1:nbParamsToUse);
    paramNames=parameterDescriptions(1:nbParamsToUse,1).';
    
    %use the directory of this program to find the .json file. 
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=fullfile(configPath,'OpenSPIM.json');
    setupConfig=loadjson(defaultConfigFileName);
    %those two values will be over written later, can be ignored. 
    setupConfig.excitation.fractionOfNumericalApertureUsed=1.0; % Full aperture can be used, depending on model
    setupConfig.detection.fractionOfNumericalApertureUsed=1.0; % Entirely open for this measurement
%     % For 10X Nikon water immersion
%     setupConfig.detection.objective.numericalAperture=0.3;
%     setupConfig.detection.objective.magnification=10;

    diffractionLimitIllumination=setupConfig.excitation.wavelength/2/(setupConfig.excitation.objective.numericalAperture*setupConfig.excitation.fractionOfNumericalApertureUsed);
    diffractionLimitDetection=setupConfig.detection.wavelength/2/(setupConfig.detection.objective.numericalAperture*setupConfig.detection.fractionOfNumericalApertureUsed);
    %real distance between each pixel based on the tube lens used. 
    dx=setupConfig.detector.pixelSize(1)/(setupConfig.detection.objective.magnification*setupConfig.detection.tubeLength/setupConfig.detection.objective.tubeLength);
    
    %read beam profile image and normalize the data, for 8bit
    %and warning for satured. 
    lightSheetMeasurement=mean(single(imread(lightSheetFileName)),3);
    lightSheetMeasurement=lightSheetMeasurement.'./(2^8-1);
    if sum(lightSheetMeasurement(:)==1)>1
        logMessage('Careful, image satured in %d pixels!',sum(lightSheetMeasurement(:)==1));
    end
    imgSize=size(lightSheetMeasurement);
    
    xRange=dx*([1:imgSize(1)]-floor(imgSize(1)/2)-1);
    zRange=round(stepSize*([1:imgSize(2)]-floor(imgSize(2)/2)-1)/1e-9)*1e-9;

%     nuZ=xRange./dx./imgSize(1)./dx.*diffractionLimitDetection;
%     OTF=(abs(nuZ)<1).*2/pi.*(acos(abs(nuZ))-abs(nuZ).*sqrt(1-nuZ.^2));
%     NSR=nuZ./3;
%     filter=OTF./(abs(OTF).^2+NSR.^2);
%     lightSheetFiltered=ifft(fft(lightSheetMeasurement).*repmat(ifftshift(filter).',[1 imgSize(2)]));
%     lightSheetFiltered=lightSheetMeasurement; % TODO remove
    
    % Remove the shear in the light sheet recording
    imageShift=zRange;
    lightSheetDesheared=lightSheetMeasurement;
    for zIdx=1:imgSize(2),
        lightSheetDesheared(:,zIdx)=interp1(xRange,lightSheetMeasurement(:,zIdx),xRange-imageShift(zIdx),'*pchip',0);
    end
    
    recordedSections=zeros(numel(xRange),numel(zSectionPositions));
    for sectionIdx=1:numel(zSectionPositions),
        recordedSections(:,sectionIdx)=mean(lightSheetDesheared(:,abs(zRange-zSectionPositions(sectionIdx))<sectionWidth/2),2);
    end
    recordedSectionsFft=fft(recordedSections-repmat(min(recordedSections),[size(recordedSections,1) 1]));
    
    %
    % Optimization
    %
    % Global optimization
    function vec=paramStructToVector(p)
        vec=zeros(1,nbParamsToUse);
        for paramIdx=1:numel(paramNames),
            vec(paramIdx)=getfield(p,paramNames{paramIdx});
        end
    end
    if globalOptimization,
        optimInfo.title = 'Airy light sheet fitting';
        objFctHandleLong=@(p)getFitError(xRange,zSectionPositions,setupConfig,recordedSectionsFft,paramStructToVector(p));
        objFctHandle=@(p)objFctHandleLong(p);
        paramDefCell=parameterDescriptions(1:nbParamsToUse,:);
        % get default DE parameters
        DEParams = getdefaultparams();
        % set number of population members (often 10*D is suggested) 
        DEParams.NP = 10*numel(params0);
        % do not use slave process here
        DEParams.feedSlaveProc = 0;
        % set times
        DEParams.maxiter       = Inf;
        DEParams.maxtime       = globalOptimizationMinutes*60;  % in seconds
        DEParams.maxclock      = [];
        % set display options
%         DEParams.refreshiter   = 1;
        DEParams.infoIterations = 1;    %new name of refreshiter
%         DEParams.refreshtime   = 10;  % in seconds
        DEParams.infoPeriod = 10;       %new name of refreshtime
%        DEParams.refreshtime2  = 20;  % in seconds
%         DEParams.refreshtime3  = 40;  % in seconds
        DEParams.sendMailPeriod = 40;   %New name of refreshtime3

        DEParams.slaveFileDir = 'C:\Windows\Temp';
        DEParams.playSound = false;
%         DEParams.saveHistory = false;
        
        % do not send E-mails
        emailParams = [];
        % set random state in order to always use the same population members here
        rand('state', 1);
        % start differential evolution
        [params, bestval]=differentialevolution(DEParams,paramDefCell,objFctHandle,{},{},emailParams,optimInfo);
    else
        params=params0;
    end
    
    params0=params;
    
    % Linear optimization
    if localOptimization,
        params=fminsearch(@(p) getFitError(xRange,zSectionPositions,setupConfig,recordedSectionsFft,p),params0,...
            optimset('Display','iter'));
    else
        params=params0;
    end
    
    %
    % Output
    %
    [err shiftAt0]=getFitError(xRange,zSectionPositions,setupConfig,recordedSectionsFft,params);
    varDescriptor=strcat('[',sprintf('%s,',paramNames{:})); varDescriptor=strcat(varDescriptor(1:end-1),'] = [');
    varDescriptor=strcat(varDescriptor,repmat('%0.3f,',[1 numel(paramNames)])); varDescriptor=strcat(varDescriptor(1:end-1),']');
    logMessage(strcat('Error of %0.6f%% for parameters ',varDescriptor,' and x-offset = %0.3f um'),[err*100 params(:).' shiftAt0*1e6]);
    fittedSections=getFittedSections(xRange,zSectionPositions,setupConfig,recordedSectionsFft,params);
    if globalOptimization || localOptimization
        save(outputFileName,'err','params','shiftAt0','stepSize','zSectionPositions','sectionWidth','xRange','zRange','lightSheetDesheared');
    end
    
    %
    % Display
    %
    figure;
    ax(1)=subplot(2,1,1);
    showImage(lightSheetDesheared.',-1,xRange*1e6,zRange*1e6);
    ax(2)=subplot(2,1,2);
    plot(repmat(xRange.'*1e6,[1 size(recordedSections,2)]),recordedSections);
    hold on;
    plot(repmat(xRange.'*1e6,[1 size(recordedSections,2)]),fittedSections,'--');
    set(ax(2:end),'XLim',xRange([1 end])*1e6,'YLim',[0 1]);
%     set(ax(3),'YLim',[0 max(fittedSections(:))]);
end

function [err shiftAt0]=getFitError(xRange,zSectionPositions,setupConfig,recordedSectionsFft,params)
    normalizedSections=calcNormalizedSections(xRange,zSectionPositions,setupConfig,params);
    [magnification, shiftAt0, tilt, err]=fitNormalizedSections(xRange,zSectionPositions,normalizedSections,recordedSectionsFft);
end
function [fittedSections shiftAt0]=getFittedSections(xRange,zSectionPositions,setupConfig,recordedSectionsFft,params)
    normalizedSections=calcNormalizedSections(xRange,zSectionPositions,setupConfig,params);
    [magnification, shiftAt0, tilt, err, fittedSections]=fitNormalizedSections(xRange,zSectionPositions,normalizedSections,recordedSectionsFft);
end
function normalizedSections=calcNormalizedSections(xRange,zSectionPositions,setupConfig,params)
    % Calculate the sections
    setupConfig.modulation=@(U,V) pupilFunctionModelForOpenSPIM(params,V);
    excitation=setupConfig.excitation;
    excitation.fractionOfNumericalApertureUsed=params(2);
    normalizedSections=calcLightSheetPsf([],zSectionPositions-params(3)*1e-6,xRange,0,setupConfig.excitation,setupConfig.modulation,setupConfig.sample.refractiveIndex);
    normalizedSections=permute(normalizedSections,[3 2 1]);
    % Add background
    normalizedSections=normalizedSections+params(1)./numel(xRange);
end
function [magnification, shiftAt0, tilt, err, fittedSections]=fitNormalizedSections(xRange,zSectionPositions,normalizedSections,recordedSectionsFft)
    overSamplingFactor=4;
    % Move and scale the sections
    normalizedSectionsFft=fft(normalizedSections);
    crossCorrelationFft=recordedSectionsFft.*conj(normalizedSectionsFft);
    if overSamplingFactor>1,
        centerPoint=1+floor(size(crossCorrelationFft,1)/2);
        crossCorrelationFft=circshift(crossCorrelationFft,[centerPoint 0]);
        crossCorrelationFft(end*overSamplingFactor,1)=0;
        crossCorrelationFft=circshift(crossCorrelationFft,-[centerPoint 0]);
    end
    [maxValue,maxI]=max(ifft(crossCorrelationFft,'symmetric')); % assume real-valued inputs
    magnifications=1./(maxValue./mean(abs(recordedSectionsFft).^2))./overSamplingFactor;
    magnification=mean(magnifications);
    % Find the shift without wrapping around the fft domain
    detectedShifts=mod(maxI+floor(size(normalizedSections,1)/2),size(normalizedSections,1))-floor(size(normalizedSections,1)/2);
    detectedShifts=detectedShifts./overSamplingFactor; % in non-oversampled pixels
    shift=mean(detectedShifts); % in non-oversampled pixels
    if numel(zSectionPositions)>1
        tilt=sum((detectedShifts-shift).*(zSectionPositions-mean(zSectionPositions)))/sum(abs(zSectionPositions-mean(zSectionPositions)).^2);
    else
        tilt=0;
    end
    dx=diff(xRange(1:2));
    positionCorrections=-shift-tilt*(zSectionPositions-mean(zSectionPositions)); % in non-oversampled pixels
    
    shiftAt0=shift-tilt*mean(zSectionPositions); % in non-oversampled pixels
    shiftAt0=dx*shiftAt0;
    if nargout>3,
        fittedSectionsFft=repmat(magnifications,[numel(xRange) 1]).*normalizedSectionsFft.*(exp(2i*pi*((ifftshift(xRange).'./dx./numel(xRange))*positionCorrections)));
        err2=mean(abs(fittedSectionsFft-recordedSectionsFft).^2)./(mean(abs(fittedSectionsFft).^2).*mean(abs(recordedSectionsFft).^2));
        err=sqrt(mean(err2(:)));
    end
    if nargout>4,
        fittedSections=real(ifft(fittedSectionsFft));
    end
end
