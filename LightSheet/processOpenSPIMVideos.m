%
%
%
 function processOpenSPIMVideos(folderNames)
 if nargin<1
    %folderNames={'D:\Light Sheet Microscopy\LabView\recorded\2014-05-28_15-57-31_FluoBeads 0.2um 49mm'};
    %folderNames={'Z:\RESULTS\OpenSPIM\2014-05-28_15-57-31_FluoBeads 0.2um 49mm'};
    %folderNames={'Z:\RESULTS\OpenSPIM\2014-05-28_17-15-33_FluoBeads 0.2um 49mm wider iris'};
    %folderNames={'C:\Users\Tom\Dropbox\OpenSPIM\data\2014-6-13 New layout 40X\2014-06-13_12-43-41_bead stepSize185nm -2.png'};
    %folderNames={'C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-12 40X\beads stepSize185nm.avi'};
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\2014-07-03\2014-07-03_12-28-45_6 beads stepSize185nm.avi';
    folderNames='F:\2015-12-07 Decon\2.6mm\2015-12-10_10-45-39 beads min-iris 2.6 FOV stepSize325nm.avi';
    %folderNames={'C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_11-56-18_2 beads stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_12-02-31_3 beads stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_12-21-08_4 beads stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_12-23-15_5 beads stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_12-28-45_6 beads stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-04-19_9 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-23-32_10 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-28-01_11 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-42-17_12 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-49-32_13 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-03_17-56-05_14 fish  stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_10-56-04_fish stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_11-37-57_larva stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_11-48-09_larva stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_11-54-08_larva stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_11-57-24_larva stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_12-09-50_larva 0.27mm stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-04_12-25-37_larva 1.8mm stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_13-11-14_larva  01 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_13-17-15_larva  02 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_13-22-22_larva  03 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_13-34-56_larva  04 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_13-56-15_larva  05 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Documents\MATLAB\Experiment data\Zebrafish chunck and larva\2014-07-08_14-06-13_larva  06 stepSize185nm.avi';
    %folderNames='C:\Users\Zhengyi\Dropbox\OpenSPIM\data\2014-06-17 even illumination 40X\2014-06-17_12-26-33_even illumination 48.5mm -2 stepSize185nm.avi';
    %folderNames='C:\Users\Tom\Dropbox\OpenSPIM\data\2014-06-17 even illumination 40X\2014-06-17_12-32-44_even illumination 48.5mm -2 stepSize185nm.avi';
    %folderNames='Z:\RESULTS\OpenSPIM\Zebrafish\2014-06-17_20-46-00_even illumination 48.5mm -2 stepSize185nm Zfish.avi';
    %folderNames='Z:\RESULTS\OpenSPIM\Zebrafish\2014-06-17_20-59-01_even illumination 48.5mm -2 stepSize185nm Zfish.avi';
    %folderNames='Z:\RESULTS\OpenSPIM\Zebrafish\2014-06-18_14-28-11_Zebrafish stepSize185nm reduced iris.avi';
    %folderNames='Z:\RESULTS\OpenSPIM\Zebrafish\2014-06-18_14-47-28_Zebrafish stepSize185nm reduced iris.avi';
 end
if ~iscell(folderNames)
    folderNames={folderNames};
end

reprocess=true;

%Use the path of current m file to find the json file under the same folder
functionName=mfilename();
configPath=mfilename('fullpath');
configPath=configPath(1:end-length(functionName));
defaultConfigFileName=fullfile(configPath,'openSPIM.json');
setupConfig=loadjson(defaultConfigFileName);
%Add some specific settings
setupConfig.detection.fractionOfNumericalApertureUsed=0.3;


for folderIdx=1:numel(folderNames),
    folderName=folderNames{folderIdx};
    [~,~,ext]=fileparts(folderName);
    
    if isempty(ext) || ~strcmp(ext,'.')
        % Process this folder unless it starts with a period
        processOpenSPIMVideo(folderName,setupConfig,reprocess);
        % Recurse through sub-folders
        folderNames=dir(folderName);
        folderNames=folderNames([folderNames.isdir]);
        folderNames=cellfun(@(e)fullfile(folderName,e),{folderNames.name},'UniformOutput',false);
        processOpenSPIMVideos(folderNames);
    end
end
end


function processOpenSPIMVideo(fileName,setupConfig,reprocess)
if exist(fileName,'file') % || exist(fullfile(fileName,'0.png'),'file')
    [pathstr,name,ext]=fileparts(fileName);
    if ~exist(strcat(pathstr,'\',name),'dir')
        mkdir(strcat(pathstr,'\',name));
    end
    outputFileName=strcat(pathstr,'\',name,'\',name,'.mat');
    %     outputFileName=strcat(fileName(1:end-4),'.mat');
    if reprocess || ~exist(outputFileName,'file')
        experimentConfigFileName=strcat(fileName(1:end-4),'.json');
        %if json file for this video exists, load,
        %or use defalt configuration (openSPIM.json)
        if (exist(experimentConfigFileName,'file'))
            setupConfig=loadjson(experimentConfigFileName);
        else
%            setupConfig=struct();
%            setupConfig.detector=struct();
            setupConfig.detector.center=[0 0]*1e-6;
            setupConfig.detector.scanShear=[.0265 .0344]; % [0  0.09]; % [fraction] vertical horizontal
            setupConfig.detector.perspectiveScaling=[659.7394 718.1654]; % [0  500]; % [m^-1] vertical horizontal
%            setupConfig.detector.scanShear=[0 0]; % [0  0.09]; % [fraction] vertical horizontal
%           setupConfig.detector.perspectiveScaling=[0 0]; % [0  500]; % [m^-1] vertical horizontal
%            setupConfig.modulation=struct();
            %experimentConfig.modulation.params=[0.155,0.912,-100.062,5.929,6.110,-6.092,0.276,-0.931,-0.932]; % and x-offset = 4.035 um
            setupConfig.modulation.params=[0.1,0.90,0,8.75,0,0,0,0,0];
            setupConfig.excitation.fractionOfNumericalApertureUsed=setupConfig.modulation.params(2);
            setupConfig.detector.center=-[0 0+0]*1e-6;
%            setupConfig.excitation=struct();
        end
            
        savejson([],setupConfig,experimentConfigFileName);
%          setupConfig=structUnion(defaultConfig,experimentConfig);
        
        if ~isa(setupConfig.modulation,'function_handle')
            setupConfig.modulation=@(U,V) pupilFunctionModelForOpenSPIM(setupConfig.modulation.params,V);
        end
        
        logMessage('Processing file %s...',fileName);
        recordedImageStack=readDataCubeFromAviFile(fileName);
        recordedImageStack=medianFilter(recordedImageStack);
        stepSize=regexp(fileName,'stepSize\s?([\d\.]+)\s?nm','tokens');
        if ~isempty(stepSize)
            stepSize=str2double(stepSize{end})*1e-9;
        else
            stepSize=1.0e-6;
            logMessage('Could not detect stepSize from file name, defaulting to %0.3f micron.',stepSize*1e6);
        end
        setupConfig.stagePositions.target=([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2)).'*[0 0 stepSize];
        % Deconvolve
        logMessage('Starting image reconstruction...');
        [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,setupConfig);
        save(outputFileName,'recordedImageStack','restoredDataCube','xRange','yRange','zRange','tRange','ZOtf','lightSheetPsf','lightSheetOtf','lightSheetDeconvFilter');
        
        %
        % Display
        %
        fig=figure('Position',[100 100 600 800]);
        projectionRawImg=squeeze(max(recordedImageStack,[],1)).';
        projectionDeconvolvedImg=squeeze(max(restoredDataCube,[],1)).';
        axs(1)=subplot(2,1,1);
        showImage(projectionRawImg,-1,yRange*1e6,zRange*1e6);
        axis image;
        labs(1,1)=xlabel('x (\mum)'); labs(1,2)=ylabel('z (\mum)');
        axs(2)=subplot(2,1,2);
        showImage(projectionDeconvolvedImg,-1,yRange*1e6,zRange*1e6);
        axis image;
        labs(2,1)=xlabel('x (\mum)'); labs(2,2)=ylabel('z (\mum)');
        set(axs(:),'LineWidth',2,'FontWeight','bold','FontSize',16,'TickDir','out','Box','off','YLim',[-1 1]*35);
        set(labs(:),'FontWeight','bold','FontSize',18);
        imwrite(projectionRawImg./max(projectionRawImg(:)),strcat(pathstr,'\',name,'\',name,'_raw.png'));
        imwrite(projectionDeconvolvedImg./max(projectionDeconvolvedImg(:)),strcat(pathstr,'\',name,'\',name,'_deconvolved.png'));
        plot2svg(strcat(outputFileName(1:end-4),'_raw_deconvolved.svg'),fig,'png');
        
    else
        logMessage('Already created %s!',outputFileName);
    end
else
    logMessage('File %s does not contain an image stack.',fileName);
end
end

% for each slice, calculate the median with 4 adjacent pixels
function dataCube=medianFilter(dataCube)
shifts=[0 0; -1 0; 1 0; 0 -1; 0 1];
for zIdx=1:size(dataCube,3),
    dataCubeSlice=dataCube(:,:,zIdx);
    for shiftIdx=1:size(shifts,1),
        dataCubeSlices(:,:,shiftIdx)=circshift(dataCubeSlice,shifts(shiftIdx,:));
    end
    dataCube(:,:,zIdx)=median(dataCubeSlices,3);
end

%     dataCube=reshape(dataCube,dataSize);
end