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

pixels = 1024;
% ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
ROIPosition = [0 (2048-pixels)/2 2048 pixels];
shortExposure = 0.01*pixels/2048;
scanningFrequency = 20;
nrStacks = 2;               % numbers of image stacks to take, with the following
                            % stackInterval to determine how long in total.
stackInterval = 1;          % Unit minutes
lateralPosition = [0];
trapDistance = 0.520;     % distance between each trap, unit mm.
folderName = 'C:\zebrafish\control';


%% initialisation
% %% laser initialization
% laser = finesse();

%% stage initialization
if 0
stage = instrfind('Type', 'serial', 'Port', 'COM12', 'Tag', '');
if isempty(stage)
    stage = serial('COM12');
else
    fclose(stage);
    stage = stage(1);
end
stage.BaudRate = 19200;
stage.Parity='none';
stage.Terminator='CR/LF';
fopen(stage);
fprintf(stage, '1VA0.1');% I suppose the speed cannot be too fast, in case sample drops
end
% fprintf(stage, '1PA4');

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
src.ExposureTime = shortExposure;
vid.TriggerRepeat=nrStacks*length(lateralPosition);
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
% vid.FramesPerTrigger = 200;
vid.FramesPerTrigger = 1/shortExposure/scanningFrequency;

% each trigger will include 2 full cycle, i.e. four volumes of images
nrCycles = 1;

frames = uint16(zeros(pixels,2048,vid.FramesPerTrigger));

absPosition = query(stage,'1TP');
absPosition = str2double(absPosition(1:7));


start(vid);


for i = 1:nrStacks
    
    tic;
    board.off;

    for n = 1:length(lateralPosition)
        newPosition = absPosition + lateralPosition(n)*trapDistance;
        fprintf(stage,strcat('1PA',num2str(newPosition)));
        pause(15);
%         laser.open;
board.laser_on;
%         start(vid);
        pause(0.1);
        trigger(vid);
        pause(2);
board.laser_off;
frames = getdata(vid);
        frames = squeeze(frames);
        c = clock;
        
        %save data into images
        fileName = strcat(folderName,'\stacks\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6))));
        saveMatrixData2ImageStack(frames,fileName);
        %show the projection and save the projection;
        stackProjection = double(max(frames,[],3));
        imshow(stackProjection/max(stackProjection(:)));
        imwrite(uint8(stackProjection),strcat(fileName,'.tiff'),'tiff');
%         stop(vid);
    end
    board.on;
    fprintf(stage,strcat('1PA',num2str(absPosition)));

    elapsedTime = toc
    %
    %     nrFrames = round(stackInterval * 60 - elapsedTime) * scanningFrequency ;
    %     src.ExposureTime = 1/scanningFrequency;
    %     vid.TriggerRepeat=inf;
    %     triggerconfig(vid,'immediate');
    % 	vid.FramesPerTrigger=1;
    %     start(vid);
    %     for n = 1:nrFrames
    %         img = getsnapshot(vid);
    %         c = clock;
    %         fileName = strcat(folderName,'\',num2str(c(4)),'.',num2str(c(5)),'.',num2str(round(c(6)*100)),'.',num2str(1000+n),'.tiff');
    %         imwrite(img,fileName,'tiff');
    %     end
    
    pause(round(stackInterval * 60 - elapsedTime));
    
end

% laser.clear;
board.clear;
fclose(stage);
delete(vid);
%
