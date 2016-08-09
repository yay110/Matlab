% This program only controls the camera to acquire images when both the
% tunable lens and the scanning mirror are running and calibrated
% The triggering cable of the camera needs to be connected to the tunable
% lens.

% This programme is for long-term imaging, it takes projection images
% continuously (same frame rate to the scanning rate), then every minutes,
% it will take a 3D image stack, do maximum intensity projection and save. 

% Created by Zhengyi Yang (zy6@st-andrews.ac.uk) on 08/08/2016

pixels = 512;
shortExposure = 0.005;
scanningFrequency = 2;
nrStacks = 2;             % numbers of image stacks to take. 
stackInterval = 0.1;      % Unit minutes
folderName = 'E:\longterm imaging testing';

%% initialisation

% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');

src = getselectedsource(vid);

vid.ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
%src.ExposureTime = 0.01/2048*pixels; %in seconds

% each trigger will include 2 full cycle, i.e. four volumes of images
nrCycles = 1;


frames = uint16(zeros(pixels,pixels,vid.FramesPerTrigger));
%newFrames = frames;

for i = 1:nrStacks
    tic;
    triggerconfig(vid,'manual');
    src.ExposureTime = shortExposure;
    vid.TriggerRepeat=0;
    vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
    start(vid);
    trigger(vid);
    pause(2);
    frames = getdata(vid);
    frames = squeeze(frames);
    c = clock;
    fileName = strcat(folderName,'\stacks\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6))));
    saveMatrixData2ImageStack(frames,fileName);
    stackProjection = double(max(frames,[],3));
    imshow(stackProjection/max(stackProjection(:)));
    imwrite(uint8(stackProjection),strcat(fileName,'.tiff'),'tiff');
    elapsedTime = toc;
    
    nrFrames = round(stackInterval * 60 - elapsedTime) * scanningFrequency ;
    src.ExposureTime = 1/scanningFrequency;
    vid.TriggerRepeat=inf;
    triggerconfig(vid,'immediate');
	vid.FramesPerTrigger=1;
    start(vid);
    for n = 1:nrFrames
        img = getsnapshot(vid);
        c = clock;
        fileName = strcat(folderName,'\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6)*100)),'.',num2str(1000+n),'.tiff');
        imwrite(img,fileName,'tiff');
    end   
    stop(vid);
end
%% To acquire images and display

% 
% for i=1:nrFigures
%     trigger(vid);
%     pause(2);
%     frames = getdata(vid);
%     frames = squeeze(frames);
%     c = clock;
%     fileName = strcat('E:\Ciona 20 hours\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6))));
%     saveMatrixData2ImageStack(frames,fileName);
%     pause(45);
% end

delete(vid);
%
