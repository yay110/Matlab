% This program only controls the camera to acquire images when both the
% tunable lens and the scanning mirror are already running and calibrated

% This programme is for long-term imaging, it could takes projection images
% continuously (same frame rate to the scanning rate), then every minutes,
% it will take a 3D image stack, do maximum intensity projection and save.

% Because there are multiple wells in the acoustic trap, when there are
% multiple objects trapped, define the relative well in laterPostions, the
% connected acuator will traslate the sample between those postions and
% does high-througput imaging.

% Created by Zhengyi Yang (zy6@st-andrews.ac.uk) on 08/08/2016

pixels = 2048;
scanningFrequency = 20;
% ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];



%define number of stacks to take and
nrStacks = 120;               % numbers of image stacks to take, with the following
% stackInterval to determine how long in total.
stackInterval = 1;          % Unit minutes
stackLength = 0.15;           %Unit minutes
folderName = 'D:\2017-03-29 Zebrafish heartbeat\3rd sample\2. adding the  drug and after';


%% pump / Arduino board initialization
% Find a serial port object.
board = Arduino('COM17');
board.connect;
board.pump_on;
% %examples about pump control
% pump.on;
% pump.off;


%% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
src = getselectedsource(vid);
triggerconfig(vid,'manual');
vid.ROIPosition = ROIPosition;
%src.ExposureTime = 0.01/2048*pixels; %in seconds
src.ExposureTime = 1/scanningFrequency;
vid.TriggerRepeat=nrStacks;
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
vid.FramesPerTrigger = stackLength*60*scanningFrequency;

frames = uint16(zeros(pixels,2048,vid.FramesPerTrigger));


start(vid);

for i = 1:nrStacks
    
    tic;
    board.pump_off;
    
    board.laser_on;
    pause(7);
    
    trigger(vid);
    pause(stackLength*60);
    board.laser_off;
    board.pump_on;
    frames = getdata(vid);
    frames = squeeze(frames);
    c = clock;
    
    %save data into images
    fileName = strcat(folderName,'\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6))));
    saveMatrixData2ImageStack(frames,fileName);
    
    stackProjection = double(max(frames,[],3));
    imshow(stackProjection/max(stackProjection(:)));
    imwrite(uint8(stackProjection),strcat(fileName,'.tiff'),'tiff');    %show the projection and save the projection;
    
    i
    elapsedTime = toc
 
    
    pause(round(stackInterval * 60 - elapsedTime));
    
end

%
% for i = 1:nrStacks
%
%     tic;
%     board.pump_off;
%
%     board.laser_on;
%     pause(7);
%
%     trigger(vid);
%     pause(stackLength*60);
%     board.laser_off;
%     frames = getdata(vid);
%     frames = squeeze(frames);
%     c = clock;
%
%     %save data into images
%     fileName = strcat(folderName,'\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6))),'pump off');
%     saveMatrixData2ImageStack(frames,fileName);
%     %show the projection and save the projection;
%
%     i
%     elapsedTime = toc
%     %
%     %     nrFrames = round(stackInterval * 60 - elapsedTime) * scanningFrequency ;
%     %     src.ExposureTime = 1/scanningFrequency;
%     %     vid.TriggerRepeat=inf;
%     %     triggerconfig(vid,'immediate');
%     % 	vid.FramesPerTrigger=1;
%     %     start(vid);
%     %     for n = 1:nrFrames
%     %         img = getsnapshot(vid);
%     %         c = clock;
%     %         fileName = strcat(folderName,'\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6)*100)),'.',num2str(1000+n),'.tiff');
%     %         imwrite(img,fileName,'tiff');
%     %     end
%
%     pause(round(stackInterval * 60 - elapsedTime));
%
% end


board.pump_off;
% laser.clear;
board.disconnect;
delete(vid);
%
