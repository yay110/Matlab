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
clear beatFrequency;

pixels = 1024;
scanningFrequency = 20;
% ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];



%define number of stacks to take and
nrStacks =52;               % numbers of image stacks to take, with the following
% stackInterval to determine how long in total.
stackInterval = 1;          % Unit minutes
stackLength = 0.25;           %Unit minutes
folderName = 'F:\2017-05-11 Verapaamil\5th sample\3. washing away drug';


%% pump / Arduino board initialization
% Find a serial port object.
board = Arduino('COM7');
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
frameNumber = scanningFrequency * stackLength * 60;

h = figure;
set(h, 'Position', [50 300 1800 500])

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
    subplot(1,3,1)
    imshow(stackProjection/max(stackProjection(:)));
    
    firstFrame = squeeze(frames(:,:,1));
    [x y] = find(firstFrame ==max(firstFrame(:)));
    
    for m = 1:frameNumber
        currentFrame = squeeze((frames(:,:,m)));
        profile(m) = sum(sum(currentFrame(x-50:x,y-50:y)));
    end
    
    subplot(1,3,2)
    plot(profile);
    
    L = length(profile);       % Length of signal
    
    Y = fft(profile);
    P2 = abs(Y/L);
    P1(:) = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = scanningFrequency*(0:(L/2))/L;
    peaks = findpeaks(P1(10:end));
    beatFrequency(i) = f(find(P1(:)==max(peaks)))
    normP1 = P1/max(peaks);
    
    subplot(2,3,3)
    plot(f, normP1);
    title(num2str(beatFrequency(i)));
    xlim([0 6]);
     ylim([0 1]);
     
     subplot(2,3,6)
     plot(beatFrequency);
     xlim([0 nrStacks]);
     ylim([0 4])
    
    saveas(h,strcat(fileName,'.tiff'),'tiff');    %show the projection and save the projection;
    
    
    i
    elapsedTime(i) = toc;
    
    pause(round(stackInterval * 60 - elapsedTime));
    
end

close all;
board.pump_off;
% laser.clear;
board.disconnect;
delete(vid);
%
