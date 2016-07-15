% This program only controls the camera to acquire images when both the
% tunable lens and the scanning mirror are running and calibrated
% The triggering cable of the camera needs to be connected to the tunable
% lens. 

% This program can either (defined by 'interpolated') interpolate the axial
% positon of the stack based on sin function (slow process) or only display
% the projection of the raw data (for alignment purpose);

% Created by Zhengyi Yang (zy6@st-andrews.ac.uk) on 15/07/2016

pixels = 256;
signalFrequency = 5;
interpolated =1;

%% initialisation

% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');

triggerconfig(vid,'hardware','risingEdge','EdgeTrigger');
src = getselectedsource(vid);

vid.ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
src.ExposureTime = 0.01/2048*pixels; %in seconds

% each trigger will include 2 full cycle, i.e. four volumes of images
nrCycles = 2;
vid.FramesPerTrigger=1/signalFrequency/src.ExposureTime*nrCycles;
vid.TriggerRepeat=0;


frames = uint16(zeros(pixels,pixels,vid.FramesPerTrigger));
newFrames = frames;

disp('Camera Ready');


%% To acquire images and display

offsetFrame = 9 ;
xRange = (1:pixels);
yRange = (1:pixels);
zRange = cumsum(abs([0,diff(sin(((1-offsetFrame):(vid.FramesPerTrigger-offsetFrame))*(2*nrCycles*pi/vid.FramesPerTrigger)))]));
zRangeOrigin = zRange(end)/vid.FramesPerTrigger:zRange(end)/vid.FramesPerTrigger:zRange(end);

start(vid);
for i=1:1
    %    trigger(vid);
    frames = getdata(vid);
    frames = squeeze(frames);
    
    % Since the axial scannning position is not uniform but follow a sin
    % wave pattern, the axial positon needs to be interpolated, based on
    % the offsetFrame number which is different for each trigger
    % Probably be able to be calibrated.
    if interpolated
        for x = 1:pixels
            for y =1:pixels
                newFrames(x,y,:) = interp1(zRange,squeeze(double(frames(x,y,:))),zRangeOrigin,'pchip');
            end
        end
        
    subplot(4,4,[1,2]);imagesc(zRange,yRange,squeeze(256*max(frames,[],1)/max(max(max(frames,[],1)))));%axis image;
    title('YZ projection of raw image');
    subplot(4,4,[5,6]);imagesc(zRangeOrigin,yRange,squeeze(256*max(newFrames,[],1)/max(max(max(newFrames,[],1)))));%axis image;
    title('YZ projection after interploation');
    subplot(4,4,[9,10]);imagesc(zRange,xRange,squeeze(256*max(frames,[],2)/max(max(max(frames,[],2)))));
    title('XZ projection of raw image');
    subplot(4,4,[13,14]);imagesc(zRangeOrigin,xRange,squeeze(256*max(newFrames,[],2)/max(max(max(newFrames,[],2)))));
    title('XZ projection after interploation');
    subplot(4,4,[3,4,7,8,11,12,15,16]);imagesc(yRange,xRange,squeeze(256*max(frames,[],3)/max(max(max(frames,[],3)))));axis image;
    title('XY projection');
    else
        subplot(2,2,1);imagesc(squeeze(256*max(frames,[],1)/max(max(max(frames,[],1)))));%axis image;
        title('YZ projection of raw image');
        subplot(2,2,3);imagesc(squeeze(256*max(frames,[],2)/max(max(max(frames,[],2)))));%axis image;
        title('XZ projection of raw image');
        subplot(2,2,[2,4]);imagesc(squeeze(256*max(frames,[],3)/max(max(max(frames,[],3)))));axis image;
        title('XY projection of raw image');
%         h(1)=subplot(2,2,1);imagesc(squeeze(256*max(frames,[],1)/max(max(max(frames,[],1)))));
%         h(2)=subplot(2,2,2);imagesc(squeeze(256*max(newFrames,[],1)/max(max(max(newFrames,[],1)))));
%         linkaxes(h,'x','y');
    end
end

delete(vid);

fileName = 'C:\acquisitionData\tif';
saveMatrixData2ImageStack(newFrames,fileName);