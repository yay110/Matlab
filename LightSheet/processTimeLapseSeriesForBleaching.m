%
% processTimeLapseSeriesForBleaching(folderName)
%
function processTimeLapseSeriesForBleaching(folderName)
    if nargin<1 || isempty(folderName)
%         folderName='H:\RESULTS\TEST\2013-12-13 13_26_15.044_withCapOnCamera';
%         folderName='H:\RESULTS\TEST\2013-12-13 14_52_18.954';
%         folderName='A:\RESULTS\TimeLapse\2013-12-13 14_58_14.851';
%         folderName='A:\RESULTS\TimeLapse\2013-12-13 16_01_41.937';
%         folderName='A:\RESULTS\TimeLapse\2013-12-15 12_47_08.601_700nmRedBeads_ranDryAtEnd'; % 200 Bessel5 beads, ran dry at the end
%         folderName='A:\RESULTS\TimeLapse\2013-12-15 14_35_13.591_700nmRedBeads'; % 100 Gauss
%         folderName='A:\RESULTS\TimeLapse\2013-12-15 17_15_37.894_rhodamine';
%         folderName='A:\RESULTS\TimeLapse\2013-12-15 18_38_28.266_rhodamine';
%         folderName='A:\RESULTS\TimeLapse\2013-12-15 18_57_54.539_fluorescein_particlesOnPDMS';
%         folderName='A:\RESULTS\TimeLapse\2013-12-16 11_06_41.497';
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 13_29_54.566_Gaussian'; % Gaussian
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 13_51_43.828_Airy'; % Airy
        
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 15_04_52.016_Bessel5'; % Bessel5
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 15_29_13.222_Bessel10';
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 16_31_34.007_Airy';
%         folderName='A:\RESULTS\TimeLapse\MCF10A_1\2013-12-16 16_50_46.088_Gaussian';

%         folderName='A:\RESULTS\TimeLapse\2013-12-16 18_57_24.858_Gaussian_WGARed';
%         folderName='A:\RESULTS\TimeLapse\2013-12-16 19_37_32.514_Bessel5_WGARed';

%         folderName='A:\RESULTS\TEST\2013-12-17 12_08_45.851';

%         folderName='A:\RESULTS\TimeLapse\Cy3Acry1MicronScans\2013-12-18 11_10_58.560_Gaussian';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry1MicronScans\2013-12-18 11_14_19.153_Airy';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry1MicronScans\2013-12-18 11_18_13.745_B10';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry1MicronScans\2013-12-18 11_22_57.206_B5';

%         folderName='A:\RESULTS\TimeLapse\Cy3Acry3MicronScan\2013-12-18 12_12_05.201_B5';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry3MicronScan\2013-12-18 12_17_13.418_B10';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry3MicronScan\2013-12-18 12_21_46.911_Airy'; % did not move position in sample
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry3MicronScan\2013-12-18 12_32_55.156_Airy'; 
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry3MicronScan\2013-12-18 12_26_34.940_Gaussian';

%         folderName='A\RESULTS\TimeLapse\Cy3Acry2MicronScan\2013-12-18 12_39_06.163_Gaussian';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry2MicronScan\2013-12-18 12_42_49.356_Airy';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry2MicronScan\2013-12-18 12_51_55.295_B5';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry2MicronScan\2013-12-18 12_46_53.885_B10_dirt';
%         folderName='A:\RESULTS\TimeLapse\Cy3Acry2MicronScan\2013-12-18 13_54_32.596_B10';

        folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide1\2013-12-19 14_36_43.331_Gaussian';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide1\2013-12-19 14_29_55.706_Airy';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide1\2013-12-19 14_21_09.446_B10';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide1\2013-12-19 14_13_02.983_B5';

%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_45_57.418_Airy1DCyl';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_31_54.440_Airy1DCyl_powerTooLow';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_25_04.860_Airy1D';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_17_04.333_AirySameAODAsGaussian';

        % Good set of measurements
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_10_46.083_Gaussian';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 16_04_07.682_Airy';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 15_55_27.876_B5';
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 15_48_44.101_B10';

%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 15_39_44.292_B10_nearTopSurface'; %near surface
%         folderName='A:\RESULTS\TimeLapse\dayOldAcrylamide\2013-12-19 15_32_50.091_B5_nearTopSurface'; %near surface

%         folderName='';
%         folderName='';
%         folderName='';
%         folderName='';
%         folderName='';
%         folderName='A:\RESULTS\TimeLapse\firstDayAcrylamide\2013-12-19 20_10_59.480_B5';
    end
    
    nbBinCols=32;
    centralCubeSize=[30 10 30]*1e-6; % swipe, prop, scan (axial)
        
    files=dir(fullfile(folderName,'recording*.avi'));
    
    names={files.name};
    if isempty(names)
        logMessage('No files found in %s.',folderName);
        return;
    end
    tokens=regexp(names,'recording([0-9]+)_lambda([0-9]+)nm_alpha([0-9]+)_beta([0-9]+)\.','tokens');
    measurements=struct('name',{},'sequenceNumber',0,'wavelength',0,'alpha',0,'beta',0,'startDateTime',0,'excitationPowers',[]);
    for nameIdx=1:numel(names),
        fullFileName=fullfile(folderName,names{nameIdx});
        measurements(nameIdx).name=fullFileName;
        
        measurements(nameIdx).sequenceNumber=str2double(tokens{nameIdx}{1}{1});
        measurements(nameIdx).wavelength=str2double(tokens{nameIdx}{1}{2});
%         measurements(nameIdx).alpha=str2double(tokens{nameIdx}{1}{3});
%         measurements(nameIdx).beta=str2double(tokens{nameIdx}{1}{4});
    end
    
    [~,lastNonZeroPowerIdx]=sort([measurements(:).sequenceNumber]);
    measurements=measurements(lastNonZeroPowerIdx);
    
    % Load data
    for nameIdx=1:numel(names),
        name=measurements(nameIdx).name;
        logMessage('Processing %s...',name);
        [dataCube maxValue]=readDataCubeFromFile(name);
        dataSize=size(dataCube);
        
        if nameIdx==1,
            integratedDataCube=zeros([1 dataSize(2) numel(names)]);
        end
        
        % Find the units
        [xRange yRange zRange measurements(nameIdx).alpha, measurements(nameIdx).beta measurements(nameIdx).startDateTime measurements(nameIdx).excitationPowers]=loadSettings(name,dataSize);
        xRange=xRange*nbBinCols;
        
        % Decide where to look
        xSel=-centralCubeSize(1)/2<=xRange & xRange<centralCubeSize(1)/2;
        zSel=-centralCubeSize(2)/2<=zRange & zRange<centralCubeSize(3)/2;
        
        if any(dataCube(:)==1)
            logMessage('Saturation at %d pixels!',sum(dataCube(:)==1));
        end
        
        integratedDataCube(1,:,nameIdx)=mean(mean(dataCube(xSel,:,zSel),3));
        
    end
    for nameIdx=1:numel(names),
        elapsedTimes(nameIdx)=etime(datevec(measurements(nameIdx).startDateTime),datevec(measurements(1).startDateTime));
        excitationPowers(nameIdx)=mean(measurements(nameIdx).excitationPowers([1 end]));
    end
    
    if any(excitationPowers(end)==0)
        zeroPower=excitationPowers==0;
        logMessage('Ignoring power measurement for %d/%d scans.',[sum(zeroPower) length(excitationPowers)]);
        
        % Correct for dropped power measurements
        prevPower=excitationPowers(1);
        for scanIdx=2:length(excitationPowers),
            if excitationPowers(scanIdx)==0,
                excitationPowers(scanIdx)=prevPower;
            else
                prevPower=excitationPowers(scanIdx);
            end
        end
    end
    
    medianExcitation=median(excitationPowers);
    excitationPowerError=std(excitationPowers)/medianExcitation;
    logMessage('Excitation power %0.0fuW +-%0.1fuW (%0.3f%%)',[medianExcitation*1e6 std(excitationPowers)*1e6 excitationPowerError*100]);
    
    excitationCorrection=excitationPowers(1)./excitationPowers;
    
    correctedIntegratedDataCube=repmat(permute(excitationCorrection(:),[3 2 1]),[1 size(integratedDataCube,2) 1]).*integratedDataCube*maxValue/nbBinCols;
        
    % Find the light sheet waist and set coordinate sytem accordingly
    scanIdx=[1:size(correctedIntegratedDataCube,3)]-1;
    fittedExpParamsFullFOV=zeros(size(correctedIntegratedDataCube,2),3);
    indexesAroundWaist=find(abs(yRange)<=30e-6);
    parfor yIdx=indexesAroundWaist,
        fittedExpParamsFullFOV(yIdx,:)=fitOffsetExponential(scanIdx,squeeze(correctedIntegratedDataCube(1,yIdx,:)).');
    end
    [~,waistI]=min(fittedExpParamsFullFOV(:,3).'); % Find worst bleaching point
    center=[0 yRange(waistI)];
    logMessage('Light sheet waist at %0.3fum.',center(2)*1e6);
%     center=[0 -3.5]*1e-6;
    yRange=yRange-center(2);
        
    % Process
    ySel=-centralCubeSize(2)/2<=yRange & yRange<centralCubeSize(2)/2;
    intensities=squeeze(mean(integratedDataCube(:,ySel,:),2)).';
    
    correctedIntensities=excitationCorrection.*intensities*maxValue/nbBinCols;
    fittedExpParams=fitOffsetExponential(scanIdx,correctedIntensities);
    fittedExp=offsetExponential(fittedExpParams,scanIdx);
    
    csvwrite_with_headers([folderName '_countsAtWaist.csv'],[scanIdx(:)+1 elapsedTimes(:) correctedIntensities(:)],{'scan number','elapsed time','intensity'});
    
    % Cut some slices through the last data cube
    projSwipe=mean(dataCube(xSel,:,:),1);
    projPropagation=mean(dataCube(:,ySel,:),2);
    projAxial=mean(dataCube(:,:,zSel),3);
        
    %
    % Output
    %
    fig=figure('Position',[50 50 800 600]);
    subplot(2,3,1);
    imagesc(yRange*1e6,xRange*1e6,squeeze(projAxial)); axis equal tight;
    xlabel('x (propagation) [\mum]'); ylabel('y (swipe) [\mum]');
    subplot(2,3,2);
    imagesc(zRange*1e6,xRange*1e6,squeeze(projPropagation)); axis equal tight;
    xlabel('z (axial) [\mum]'); ylabel('y (swipe) [\mum]');
    subplot(2,3,3);
    plot(yRange*1e6,fittedExpParamsFullFOV(:,3).');
    xlabel('y (propagation) [um]'); xlim([-1 1]*30); set(gca,'XTick',[-100:5:100]);
    subplot(2,3,4);
    imagesc(yRange*1e6,zRange*1e6,squeeze(projSwipe).'); axis equal tight;
    xlabel('x (propagation) [\mum]'); ylabel('z (axial) [\mum]');
    subplot(2,3,5);
%     plot(elapsedTimes/60,correctedIntensities);
%     xlabel('elapsed time [min.]'); ylabel('average counts [a.u.]');
    plot(scanIdx,correctedIntensities,scanIdx,fittedExp); ylim([min(correctedIntensities) max(correctedIntensities)]);
    xlabel('scan # [scans]'); ylabel('average counts [a.u.]');
    title(sprintf('%0.3f + %0.3f exp(%0.3f #)',fittedExpParams));
    
    set(fig,'NumberTitle','off','Name',folderName);
%     logMessage('%0.3f ',correctedIntensities);
end

function y=offsetExponential(params,x)
    y=params(1)+params(2)*exp(params(3).*x);
end

function fitParams=fitOffsetExponential(x,y)
    dydx=diff(y)./diff(x);
    ddydx=diff(dydx)./diff(x(1:end-1));
    offset=min(y);
    rate=min(0,median(ddydx./(0.5*(dydx(1:end-1)+dydx(2:end)))));
    factor=median((y-offset)./exp(rate*x));
    fitParams0=double([offset factor rate]);
    
    errorFunction=@(p)sum(abs(offsetExponential(p,x)-y).^2);
    
    fitParams=fminsearch(errorFunction,fitParams0,optimset('Display','off'));
    fitParams=fminsearch(errorFunction,fitParams,optimset('Display','off'));
    fitParams=fminsearch(errorFunction,fitParams,optimset('Display','off'));
    
%     figure('Name',sprintf('fit parameters: y = %0.3f %0.3f exp(%0.3f x)',fitParams));
%     plot(x,y,x,offsetExponential(fitParams,x),x,offsetExponential(fitParams,x)); ylim([min(y) max(y)]);

    fitOpts=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[max(y) max(y) 1],'StartPoint',fitParams0);
    fit(x,y,fittype(@(a,b,c,x,y) a+b*exp(-c*x)-y,'problem','n','options',fitOpts))
end

function [xRange yRange zRange alpha beta startDateTime excitationPowers]=loadSettings(specificConfigFileNames,dataCubeSize)
    %Load the configuration
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'/waterImmersion.json');
    defaultConfig=loadjson(defaultConfigFileName);
    if (~iscell(specificConfigFileNames))
        specificConfigFileNames={specificConfigFileNames};
    end
    % Find a json file
    fileIdx=1;
    while(fileIdx<=length(specificConfigFileNames) && ~exist([specificConfigFileNames{fileIdx}(1:end-4),'.json'],'file')),
        fileIdx=fileIdx+1;
    end
    if (fileIdx<=length(specificConfigFileNames))
        specificConfig=loadjson([specificConfigFileNames{fileIdx}(1:end-4),'.json']);
        config=structUnion(defaultConfig,specificConfig);
    else
        logMessage('Description file with extension .json is missing, trying to load .mat file.');
        try
            matFile=matfile(specificConfigFileNames{1},'Writable',false);
            config=matFile.setupConfig;
        catch Exc
            logMessage('Failed... assuming defaults!');
            config=defaultConfig;
        end
    end
    
    if (~isfield(config.detector,'center') || isempty(config.detector.center))
        config.detector.center=[0 0];
    end
    
    % Prepare the results
    stageTranslationStepSize=norm(median(diff(config.stagePositions.target)));
    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    xRange=-config.detector.center(1)+config.detector.pixelSize(1)*([1:dataCubeSize(1)]-floor(dataCubeSize(1)/2)-1)/realMagnification; % up/down
    yRange=-config.detector.center(2)+config.detector.pixelSize(2)*([1:dataCubeSize(2)]-floor(dataCubeSize(2)/2)-1)/realMagnification; % left/right
    zRange=stageTranslationStepSize*([1:dataCubeSize(3)]-floor(dataCubeSize(3)/2+1))/config.detector.framesPerSecond; %Translation range (along z-axis)
    
    alpha=config.modulation.alpha;
    beta=config.modulation.beta;
    
    startDateTime=datenum(config.startDateTime,'yyyy-mm-ddTHH:MM:SS.FFF');
    
    excitationPowers=config.excitation.powers;
end
    