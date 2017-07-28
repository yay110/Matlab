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
binning = 2;
scanningFrequency = 20;
% ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
ROIPosition = [(2048/binning-pixels)/2 (2048/binning-pixels)/2 pixels pixels];



%define number of stacks to take and
nrStacks = 60;               % numbers of image stacks to take, with the following
% stackInterval to determine how long in t otal.
stackInterval = 1;          % Unit minutes
stackLength = 0.15;           %Unit minutes
folderName = 'D:\2017-03-31 Zebrafish heartbeat\3rd sample\2. adding drug';


%% pump / Arduino board initialization
% Find a serial port object.
board = Arduino('COM17');
board.connect;
board.pump_on;
% %examples about pump control
% pump.on;
% pump.off;

%define parameters of the signal for FFT
Fs = scanningFrequency;                        %sampling frequency
T = 1/Fs;                       % sampling period
sectionLength = stackLength *60;
L = sectionLength*Fs;       % Length of signal
t = (0:L-1)*T;                  %Time vector

%% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_BIN2x2_1024x1024_FastMode');
src = getselectedsource(vid);
triggerconfig(vid,'manual');
vid.ROIPosition = ROIPosition;
%src.ExposureTime = 0.01/2048*pixels; %in seconds
src.ExposureTime = 1/scanningFrequency;
vid.TriggerRepeat=nrStacks;
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
vid.FramesPerTrigger = stackLength*60*scanningFrequency;

frames = uint16(zeros(pixels,2048,vid.FramesPerTrigger));

h=figure;

start(vid);

for i = 1:nrStacks
    
    i
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
    subplot(1,2,1)
    imshow(stackProjection/max(stackProjection(:)));
%     imwrite(uint8(stackProjection),strcat(fileName,'.tiff'),'tiff');    %show the projection and save the projection;
    
    %%============== find area of heart ===================
%     subDataCube = double(frames)./double(max(frames(:)));
%     for m = 1:size(subDataCube,3)
%         I = subDataCube(:,:,m);
%         
%         BW = im2bw(I, 0.2);
%         %             subplot(1,2,1),imagesc(I);
%         %             subplot(1,2,2),imagesc(BW);
%         n1               	= 200;
%         BW1               	= bwareaopen(BW,n1);
%         n2                  = 1;
%         BW2              	= imclose(BW1,strel('disk', n2));       %n2 = 1, close the structure
%         BW3              	= imfill(BW2,'holes');
%         %             imagesc(BW3)
%         %             hold on;        
%         B = bwboundaries(BW3,'noholes');
%         %             for k = 1:length(B)
%         %                 boundary = B{k};
%         %                 plot(boundary(:,2), boundary(:,1));
%         %             end
%         %             hold off;        
%         heartBoundary(i,m) = bwarea(BW3);
%         %             measurements = regionprops(BW3, 'Area');
%         %             heatBoundary(m) = [measurements.Area];
%     end
%     subProfile = heartBoundary(i,:);
%     
    %%================== plot Z profile ===================
    for m = 1:size(frames,3)
    subProfile(m) = sum(sum(frames(:,:,m)));
    end
    
    L = sectionLength*Fs;
    Y = fft(subProfile);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    peaks = findpeaks(P1);
    
%     beatFreqency(i) = f(find(P1==max(peaks)));
%     normP1 = P1/max(peaks);
    subplot(1,2,2)
    plot(f,P1)
    xlim([1 3])
    saveas(h,strcat(fileName,'figure.tiff'))

%     beatFreqency(i)
    
    elapsedTime = toc
    
    pause(round(stackInterval * 60 - elapsedTime));
    
end

% save processed.mat;

board.pump_off;
% laser.clear;
board.disconnect;
delete(vid);
%
