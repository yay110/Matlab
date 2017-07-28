function createSimulatedLightSheetMovie()
    close all;
    
    displayProgress=true;
    movieFileName='CoatOfArms';
    movieSize=[480 640];
    
    twoPhoton=false;
    
    % General constants
    NAvogadro=6.0221417930e23;
    hPlanck=6.6260695729e-34;
    cLight=299792458;
    
    framesPerSecond=25;
    
    %% Calculation parameters
    
    deflectBeamInsteadOfSampleMovement=false;
    
    % Choice of light-sheet
%     illuminationCroppingFactors=[1 1; 1 1]*.06; % obsolete
%     lightSheetName='Gaussian';
%     alpha=0.0;
%     openFractionOfRadius=1.0;
%     lightSheetName='Bessel10';
%     alpha=0.0;
%     openFractionOfRadius=0.1;
%     lightSheetName='Bessel5';
%     alpha=0.0;
%     openFractionOfRadius=0.05;
    lightSheetName='Airy';
    alpha=7.0;
    openFractionOfRadius=1.0;

    movieFileName=strcat(movieFileName,lightSheetName);
        
    % Objectives and wavelengths
    outsideCapillary=false;
    config=loadSetUpData(outsideCapillary);
    if (twoPhoton)
        config.excitation.wavelength=config.excitation.wavelength*2;
        config.excitation.power=config.excitation.power/5;
%         config.sample.fluorophore.maxDensity=(1/5)*config.sample.fluorophore.maxDensity;
        movieFileName=strcat('TwoPhoton',movieFileName);
    end
      
    % Sample grid specification
    xRange=[0]*1e-6; % up
    yRange=[-10:(1/3):100]*1e-6; % light-sheet propagation axis
    zRange=[-50:(1/3):50]*1e-6; % detection axis
    tRange=[-60:(1/3):60]*1e-6; %[-100:.333:100]*1e-6; %Scan range (along z-axis)
    xRange=single(xRange); yRange=single(yRange); zRange=single(zRange); tRange=single(tRange);
%     elementSize=diff([xRange(1:2); yRange(1:2); zRange(1:2)].');
    elementSize=diff([xRange([1 1]); yRange(1:2); zRange(1:2)].');
    
    % Sample definition
%     testImage=getTestImage('StAndrewsCoatOfArms128x156.png');
    testImage=getTestImage('StAndrewsCoatOfArms256x312.png');
    testImage=mean(testImage,3);
    testImage=testImage(1:1:end,1:1:end)./max(testImage(:));
    %Add noise around test object
    randStream=RandStream('mt19937ar', 'Seed',0);
    %testImage=repmat(testImage,[1 2]).';
     %Create the sample volume
    backgroundLevel=config.sample.backgroundLevel; % fraction of dynamic range
    backgroundNoise=0.0; % fraction of dynamic range
    fluorescenceDensity=max(0,backgroundLevel+backgroundNoise*randn(randStream,[length(xRange) length(yRange) length(zRange)]));
    fluorescenceDensity(1+floor(end/2), floor(end/2)+1+[1:size(testImage,1)]-floor(size(testImage,1)/2)-1, floor(end/2)+1+[1:size(testImage,2)]-floor(size(testImage,2)/2)-1)=permute(testImage,[1 3 2]);
    clear('testImage');
    fluorescenceDensity=fluorescenceDensity*config.sample.fluorophore.maxDensity; % fluorophores * m^-3
    fluorescenceDensity=single(fluorescenceDensity);
    
    %% Preparative Calculations
    excitationFocalLength=config.excitation.objective.tubeLength/config.excitation.objective.magnification;
    
    %The excitation beam travels along the x-axis
    absorptionProbabilityOfFluorophore=config.sample.fluorophore.extinctionCoefficient*elementSize(2);

    % Define some noise related variables
    detectionObjectiveSolidAngle=2*pi*(1-cos(asin((config.detection.objective.numericalAperture/config.detection.objective.refractiveIndex)/config.sample.refractiveIndex)));
    objectiveEfficiency=detectionObjectiveSolidAngle/(4*pi);
    overallDetectionEfficiency=config.sample.fluorophore.quantumYield*objectiveEfficiency*config.detector.quantumEfficiency;
    photoElectronsPerGrayLevel=config.detector.wellDepth/config.detector.numberOfGrayLevels;
    detectorNoise=config.detector.integrationTime*config.detector.darkCurrent+config.detector.readOutNoise; % e-
      
    % Pre-calc. detection PSF
%     detectionPsf=calcDetectionPsf(xRange,yRange,zRange,0,0,config.detection,config.sample.refractiveIndex); %Assume shift invariance of the detection PSF
    detectionPsf=calcVectorialPsf(xRange,yRange,zRange,config.detection.wavelength,@(U,V) 1,@(U,V) 1i,config.detection.objective.numericalAperture,config.sample.refractiveIndex,config.detection.objective.magnification,config.detection.objective.tubeLength); %Assume shift invariance of the detection PSF
    detectionPsfCenterIdx=[find(abs(xRange)==min(abs(xRange)),1) find(abs(yRange)==min(abs(yRange)),1)];
    detectionPsf=overallDetectionEfficiency*detectionPsf;
    detectionOtf=fft2(circshift(detectionPsf,[-detectionPsfCenterIdx 0]));
    
    %Pre-calc the light-sheet
    photonEnergy=hPlanck*cLight/config.excitation.wavelength;
    [psf psfTwoPhotonSwiped]=calcLightSheetPsf(xRange,yRange,zRange,0,config.excitation,alpha,openFractionOfRadius,config.sample.refractiveIndex);
%     psf=squeeze(psf);
%     psfMax=max(psf,[],2);
%     psfNorm=psf.*repmat(1./psfMax(:),[1 size(psf,2)]);
%     close all;
%     imagesc(psfNorm.'>exp(-2));
%     axis equal;
    
    if (twoPhoton)
        lightSheetPsf=psfTwoPhotonSwiped;
    else
        lightSheetPsf=psf;
    end
    lightSheetPsf=lightSheetPsf./(mean(lightSheetPsf(:))*size(lightSheetPsf,3));
    lightSheetPsf=lightSheetPsf*config.excitation.power*config.detector.integrationTime/photonEnergy; %Convert to photons per voxel
            
    %Prepare scan-progress figure
    if (displayProgress)
        lineWidth=2;
        scanFig=figure();
        scanFigAx(1)=subplot(2,3,1,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        scanFigAx(2)=subplot(2,3,2,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        scanFigAx(5)=subplot(2,3,3,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        scanFigAx(6)=subplot(2,3,6,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        if (~isempty(movieFileName))
            movieFig=figure('Color',[1 1 1],'Position',[0 0 movieSize(2)+1 movieSize(1)+1]);
            scanFigAx(3)=subplot(1,2,1,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
            scanFigAx(4)=subplot(1,2,2,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
            aspectRatio=.6; %length(xRange)/length(zRange);
            set(scanFigAx(3),'Units','pix','Position',[80 30 round(270*[1 1/aspectRatio])],'TickDir','out','Box','on');
            set(scanFigAx(4),'Units','pix','Position',[361 30 round(270*[1 1/aspectRatio])],'TickDir','out','Box','on');
            set(scanFigAx(3),'XTickLabelMode','manual','XTick',20*[-10:10],'XTickLabel',20*[-10:10]);
            set(scanFigAx(3),'YTickLabelMode','manual','YTick',20*[-10:10],'YTickLabel',20*[-10:10]);
            set(scanFigAx(4),'XTickLabelMode','manual','XTick',20*[-10:10],'XTickLabel',20*[-10:10]);
            set(scanFigAx(4),'YTickLabelMode','manual','YTickLabel',[]);
%             axes('Units','pix','Position',get(scanFigAx(3),'Position'),'TickDir','out','XAxisLocation','top','YAxisLocation','right');
%             axes('Units','pix','Position',get(scanFigAx(4),'Position'),'TickDir','out','XAxisLocation','top','YAxisLocation','right');
%             axis(scanFigAx(3),'equal');
%             axis(scanFigAx(4),'equal');

            try
                writerObj = VideoWriter(strcat(movieFileName,'.avi'),'Uncompressed AVI'); %'Motion JPEG AVI'); %'Uncompressed AVI');
                writerObj.FrameRate=framesPerSecond;
            catch Exc
                writerObj = videoWriter(strcat(movieFileName,'.avi'),'width',movieSize(2)/2,'height',movieSize(1)/2); %,'codec','xvid'); %,'fps',25);
            end
            
            open(writerObj);                
        else
            scanFigAx(3)=subplot(2,3,4,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
            scanFigAx(4)=subplot(2,3,5,'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        end

        detectionPsfNormalization=max(max(detectionPsf(:,:,floor(end/2)+1)));
    
        %Calculate detectionNormalization
        emission=lightSheetPsf.*(absorptionProbabilityOfFluorophore*fluorescenceDensity)*overallDetectionEfficiency;
        detectionNormalization=0.1*max(emission(:))/max(max(sum(detectionPsf,3)));
        
        lightSheetNormalization=max(lightSheetPsf(:));
        if (twoPhoton)
            lightSheetNormalization=5*lightSheetNormalization;
        end
    else
        scanFig=[];
    end
    
    try
        %% Scan the light-sheet
        recordedImageStack=zeros([length(xRange),length(yRange),length(tRange)],'single');
        for (tIdx=1:length(tRange))
            zIdx=find(abs(zRange-tRange(tIdx))==min(abs(zRange-tRange(tIdx))),1);

            lightSheet=zeros(size(fluorescenceDensity),'single');
            theta=-tRange(tIdx)/excitationFocalLength;
            if (deflectBeamInsteadOfSampleMovement)
                YI=repmat(cos(theta)*(yRange+excitationFocalLength).',[1 length(zRange)]) - repmat(sin(theta)*zRange,[length(yRange) 1]) - excitationFocalLength;
                ZI=repmat(sin(theta)*(yRange+excitationFocalLength).',[1 length(zRange)]) + repmat(cos(theta)*zRange,[length(yRange) 1]);
            else
                % No rotation, only translation but at the same rate
                YI=repmat(yRange.',[1 length(zRange)]);
                ZI=repmat(theta*(yRange+excitationFocalLength).',[1 length(zRange)]) + repmat(zRange,[length(yRange) 1]);
            end
            for (xIdx=1:length(xRange))
                lightSheet(xIdx,:,:)=ipermute(interp2(zRange,yRange.',permute(lightSheetPsf(xIdx,:,:),[2 3 1]),ZI,YI,'cubic',0),[2 3 1]);
    %             lightSheet(:,yIdx,:)=circshift(lightSheetPsf(:,yIdx,:),[0 0 zIdx-floor(size(lightSheetPsf,3)/2)]);
            end
            %Simulate excitation
            excitedFluorophoreCount=single(poisson(double(absorptionProbabilityOfFluorophore*fluorescenceDensity.*lightSheet)));

            %Pad-shift detection psf, and convolve
            zDisplacedDetectionOtf=detectionOtf;
            zDisplacedDetectionOtf(:,:,2*end)=0;
            zDisplacedDetectionOtf=circshift(zDisplacedDetectionOtf,[0 0 zIdx-floor(size(detectionOtf,3)/2)-1]);
            zDisplacedDetectionOtf=zDisplacedDetectionOtf(:,:,1:end/2);
            %Convolve excitation with detectionPsf
            emission=ifft2(fft2(excitedFluorophoreCount).*zDisplacedDetectionOtf,'symmetric');
            emission=max(0,emission); %Avoid negative values (rounding errors?)
            recordedImage=floor((poisson(sum(emission,3))+detectorNoise*randn([size(emission,1),size(emission,2)]))/photoElectronsPerGrayLevel); %Add noise to the recording
            if (any(recordedImage(:)>config.detector.numberOfGrayLevels-1))
                logMessage('%u pixels saturated',sum(recordedImage(:)>detector.numberOfGrayLevels-1));
            end
            recordedImage=min(max(recordedImage,0),config.detector.numberOfGrayLevels-1);
            recordedImageStack(:,:,tIdx)=recordedImage;

            %Display
            if (displayProgress)
                zDisplacedDetectionPsf=circshift(ifft2(zDisplacedDetectionOtf,'symmetric'),[detectionPsfCenterIdx 0]);

                lightSheetRecordingCameraView=squeeze(lightSheet(:,:,floor(end/2)+1)).'/lightSheetNormalization; % y,x
                lightSheetTopView=squeeze(lightSheet(floor(end/2)+1,:,:)).'/lightSheetNormalization; % z,x
                lightSheetSideView=squeeze(lightSheet(:,floor(end/2)+1,:)).'/lightSheetNormalization; % y,z
                emissionRecordingCameraView=squeeze(sum(emission,3)).'/detectionNormalization; % y,x
                emissionTopView=squeeze(sum(emission,1)).'/detectionNormalization; % z,x
                emissionSideView=squeeze(sum(emission,2)).'/detectionNormalization; % y,z
                axes(scanFigAx(1)); showImage(spectrumToRGB([config.excitation.wavelength, config.detection.wavelength],cat(3,lightSheetRecordingCameraView,emissionRecordingCameraView)),[],xRange*1e6,yRange*1e6,scanFigAx(1)); xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);ylabel('y [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); set(scanFigAx(1),'LineWidth',lineWidth,'TickDir','out');
                axis('equal');
                title('Recording Camera View','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
                axes(scanFigAx(2)); showImage(spectrumToRGB([config.excitation.wavelength, config.detection.wavelength],cat(3,lightSheetSideView,emissionSideView)),[],zRange*1e6,xRange*1e6,scanFigAx(2)); xlabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);ylabel('y [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); set(scanFigAx(2),'XDir','reverse','LineWidth',lineWidth,'TickDir','out');
                title('Side View','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);

                axes(scanFigAx(6)); showImage(squeeze(zDisplacedDetectionPsf(floor(end/2)+1,:,:)).'/detectionPsfNormalization,[],yRange*1e6,zRange*1e6,scanFigAx(6)); xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); set(scanFigAx(6),'LineWidth',lineWidth,'TickDir','out');
                title('Detection PSF - Top View','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
                axes(scanFigAx(5)); showImage(zDisplacedDetectionPsf(:,:,zIdx).'/detectionPsfNormalization,-1,xRange*1e6,yRange*1e6,scanFigAx(5)); xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);ylabel('y [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); set(scanFigAx(5),'LineWidth',lineWidth,'TickDir','out');
                title('Detection PSF','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);

                imageStackSlice=squeeze(recordedImageStack(floor(end/2)+1,:,:));
                notSaturated=double(imageStackSlice<config.detector.numberOfGrayLevels-1);
                axes(scanFigAx(4)); showImage(cat(3,imageStackSlice,imageStackSlice.*notSaturated,imageStackSlice.*notSaturated)/(config.detector.numberOfGrayLevels-1),[],tRange*1e6,yRange*1e6,scanFigAx(4));
                xlabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
                %ylabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
                set(scanFigAx(4),'LineWidth',lineWidth,'TickDir','out');
                axis('equal'); xlim(zRange([1 end])*1e6); ylim(yRange([1 end])*1e6);
                set(scanFigAx(4),'XTickLabelMode','manual','XTick',20*[-10:10],'XTickLabel',20*[-10:10]);
                set(scanFigAx(4),'YTickLabelMode','manual','YTickLabel',[]);
                title('recorded image','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
                axes(scanFigAx(3)); showImage(spectrumToRGB([config.excitation.wavelength, config.detection.wavelength],cat(3,lightSheetTopView.',emissionTopView.')),[],zRange*1e6,yRange*1e6,scanFigAx(3));
                xlabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
                ylabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
                set(scanFigAx(3),'LineWidth',lineWidth,'TickDir','out');
                axis('equal'); xlim(zRange([1 end])*1e6); ylim(yRange([1 end])*1e6);
                set(scanFigAx(3),'XTickLabelMode','manual','XTick',20*[-10:10],'XTickLabel',20*[-10:10]);
                set(scanFigAx(3),'YTickLabelMode','manual','YTick',20*[-10:10],'YTickLabel',20*[-10:10]);
                title('light sheet','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);

                drawnow();

                if (~isempty(movieFileName))
                    recordFrame(movieFig,writerObj);
                end
            end
            %Log progess in any case
            pctDone=round(100*(tIdx)/length(tRange));
            if (pctDone>round(100*(tIdx-1)/length(tRange)))
                logMessage('%u%% done',pctDone);
            end
        end

        %% Image reconstruction
        %[restoredDataCubeI lightSheetDeconvFilter lightSheetOtf ZOtf tRangeExtended]=reconstructLightSheetDataCube(xRange,yRange,zRange,tRange,recordedImageStack,config.excitation,config.detection,lightSheetPsf,detectionPsf,config.sample.signalLevel,backgroundLevel,deflectBeamInsteadOfSampleMovement);
        [restoredDataCubeI lightSheetDeconvFilter lightSheetOtf ZOtf xRange yRange zRange tRangeExtended]=deconvolveRecordedImageStack(recordedImageStack,config);

        %% Display results
        if (ishandle(scanFig))
            close(scanFig);
        end
        resultFig=figure();
        resultAx(1)=subplot(2,2,1);
        if (~isempty(movieFileName))
            resultAx(2)=scanFigAx(3);
        else
            resultAx(2)=subplot(2,2,2);
        end
        resultAx(3)=subplot(2,2,3);
        resultAx(4)=subplot(2,2,4);
        %The recorded image
        axes(resultAx(1)); showImage(squeeze(recordedImageStack(:,floor(end/2)+1,:)).'./max(abs(recordedImageStack(:))),[],xRange*1e6,tRange*1e6,resultAx(1));
        set(resultAx(1),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        title(sprintf('%s',lightSheetName),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
        %The restored image
        axes(resultAx(2)); showImage(squeeze(restoredDataCubeI(:,floor(end/2)+1,:))./max(abs(restoredDataCubeI(:))),[],tRangeExtended*1e6,xRange*1e6,resultAx(2));
        set(resultAx(2),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        xlabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        ylabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        set(resultAx(2),'LineWidth',lineWidth,'TickDir','out');
        axis('equal'); xlim(zRange([1 end])*1e6); ylim(xRange([1 end])*1e6);
        if (~isempty(movieFileName))
            title('deconvolved','FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
            set(scanFigAx(3),'XTickLabelMode','manual','XTick',20*[-10:10],'XTickLabel',20*[-10:10]);
            set(scanFigAx(3),'YTickLabelMode','manual','YTick',20*[-10:10],'YTickLabel',20*[-10:10]);
            drawnow();
            for (repFrame=1:25*5)
                recordFrame(movieFig,writerObj);
            end
            saveas(movieFig,strcat(movieFileName,'.eps'),'epsc2');
        else
            title(sprintf('%s restored',lightSheetName),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
        end
        %The light-sheet
        axes(resultAx(3)); showImage(squeeze(lightSheetPsf(:,floor(end/2)+1,:)).'./max(abs(lightSheetPsf(:))),[],xRange*1e6,zRange*1e6,resultAx(3));
        set(resultAx(3),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        title(sprintf('%s light-sheet',lightSheetName),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
        %The MTF
        axes(resultAx(4));
        plotSliceXPosIdx=floor(length(xRange)/2)+1;
    %     plotSliceXPosIdx=find(xRange==0);
        mtf=abs(squeeze(lightSheetOtf(plotSliceXPosIdx,floor(end/2)+1,[floor(end/2)+1:end,1]) ));
        filterAmplification=abs(squeeze(lightSheetDeconvFilter(plotSliceXPosIdx,floor(end/2)+1,[floor(end/2)+1:end,1])));
        mtfRestored=mtf.*filterAmplification;
        noiseLevel=0.01;
        noiseAmplification=noiseLevel*filterAmplification;
        ZOtfRange=-1e-6*ZOtf(floor(end/2)+1:-1:1);
        area(ZOtfRange,noiseAmplification,'LineWidth',lineWidth,'FaceColor',[.5 .5 .5],'EdgeColor','none');
        hold on;
        plot(ZOtfRange,mtf,'Color',[.8 0 0],'LineWidth',lineWidth,'LineStyle','-');
        plot(ZOtfRange,mtfRestored,'Color',[0 .75 .1],'LineWidth',2,'LineStyle','-');
        xlim([0 ZOtfRange(end)]);
        ylim([0 1.1]);
        hold off;
        set(resultAx(4),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',18,'LineWidth',lineWidth);
        xlabel('\nu_z [cycles/\mum]','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18); ylabel('MTF','Interpreter','Tex','FontWeight','bold','FontName','Arial Unicode MS','FontSize',18);
        title(sprintf('%s Deconvolution MTF',lightSheetName),'FontWeight','bold','FontName','Arial Unicode MS','FontSize',24);
    
    catch Exc
        if (~isempty(movieFileName))
            close(writerObj);
        end
        rethrow(Exc);
    end
    
    if (~isempty(movieFileName))
        close(writerObj);
    end
end

function recordFrame(movieFig,writerObj)
    set(movieFig,'Color',[1 1 1]);
    frm=getframe(movieFig);
    img=frame2im(frm);
    
    %Down sample for anti-aliasing
    if (mod(size(img,1),2)==0)
        img(end+1,1,1)=0;
    end
    if (mod(size(img,2),2)==0)
        img(1,end+1,1)=0;
    end
    img=double(img)./255; % convert from bytes
%     img=(img(1:2:end-1,:,:)+img(2:2:end,:,:))./2;
%     img=(img(:,1:2:end-1,:)+img(:,2:2:end,:))./2;
    img=(img(1:2:end-1,:,:)+2*img(2:2:end,:,:)+img(3:2:end,:,:))./4;
    img=(img(:,1:2:end-1,:)+2*img(:,2:2:end,:)+img(:,3:2:end,:))./4;
    
    if (ismethod(writerObj,'writeVideo'))
        % The Mathworks VideoWriter object
        writeVideo(writerObj,img);
    else
        addframe(writerObj,img);
    end
end

function showSlices(dataCube,timeInterval,xRange,yRange,zRange)
    ax=gca;
    
    if (nargin<2)
        timeInterval=0.5;
    end
    if (nargin<3)
        xRange=[];
    end
    if (nargin<4)
        yRange=[];
    end
    if (nargin<5)
        zRange=[];
    end
    
    for (zIdx=1:size(dataCube,3))
        if (ishandle(ax))
            showImage(dataCube(:,:,zIdx).',[],xRange,yRange,ax);
            if (~isempty(zRange))
                title(sprintf('slice: %d',zRange(zIdx)));
            end
            pause(timeInterval);
            drawnow;
        else
            break;
        end
    end
end


function config=loadSetUpData(outsideCapillary)
    %Load the default configuration
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'capillarySetup.json');
    config=loadjson(defaultConfigFileName);
    
%     % Objectives and wavelengths
%     config.excitation={};
%     config.excitation.wavelength=532e-9; % Firefly* % 395e-9; for GFP
%     config.excitation.objective={};
%     config.excitation.objective.numericalAperture=.42;
%     config.excitation.objective.magnification=20;
%     config.excitation.objective.tubeLength=200e-3; %for Mitutoyo the rear focal length is 200mm
    config.excitation.power=0.01e-3;
% 
%     detection={};
%     detection.wavelength=600e-9; % Firefly* % 509e-9; for GFP
%     detection.objective={};
%     detection.objective.numericalAperture=.40; % .28;
%     detection.objective.refractiveIndex=1.0; % air
%     detection.objective.magnification=22; %Actual magnification
%     detection.objective.tubeLength=1.1*160e-3; % Newport %200e-3; % Mitutoyo

    % Sample medium specification
    NAvogadro=6.0221417930e23;
    config.sample={};
    config.sample.fluorophore={};
    config.sample.fluorophore.maxDensity=30e-6*1e3*NAvogadro; % m^-3, 30umol/L density 
    config.sample.fluorophore.extinctionCoefficient=30000*(100/10^3)/NAvogadro; % m^2
    config.sample.fluorophore.quantumYield=0.79;
    config.sample.refractiveIndex=1.4; %1.33; %1.4; %Assuming the objective NA is in air
    config.sample.signalLevel=0.03; % for deconvolution
    config.sample.backgroundLevel=0.01; % fraction of dynamic range

%     % Detector specification
%     detector={};
%     detector.integrationTime=1e-3; % s
%     detector.quantumEfficiency=0.55;
%     detector.darkCurrent=200; % e-/s
%     detector.readOutNoise=16; % e-
%     detector.wellDepth=18e3; % e-
%     detector.numberOfGrayLevels=2^8;
%     detector.pixelSize=7.4*[1 1]*1e-6;
%     detector.framesPerSecond=1; % [s]
    
    if (outsideCapillary)
        config.sample.refractiveIndex=1.0;
    end
end