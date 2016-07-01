pixels = 256;
signalFrequency = 5;

%% initialisation

% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');

triggerconfig(vid,'hardware','risingEdge','EdgeTrigger');
%triggerconfig(vid,'manual');
src = getselectedsource(vid);

vid.ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
src.ExposureTime = 0.01/2048*pixels; %in seconds

% each trigger will include 2 full cycle, i.e. four volumes of images
nrCycles = 2;
vid.FramesPerTrigger=1/signalFrequency/src.ExposureTime*nrCylces;
%vid.TriggerRepeat=Inf;


frames = uint16(zeros(pixels,pixels,vid.FramesPerTrigger));
newFrames = frames;

disp('Camera Ready');

offsetFrame = 9 ;
xRange = (1:pixels);
yRange = (1:pixels);
zRange = cumsum(abs([0,diff(sin(((1-offsetFrame):(vid.FramesPerTrigger-offsetFrame))*(2*nrCycles*pi/vid.FramesPerTrigger)))]));
zRangeOrigin = zRange(end)/vid.FramesPerTrigger:zRange(end)/vid.FramesPerTrigger:zRange(end);

start(vid);

for i=1:1
    %    trigger(vid);
    frames = getdata(vid);
    
    figure;
    frames = squeeze(frames);
    for x = 1:pixels
        for y =1:pixels
            newFrames(x,y,:) = interp1(zRange,squeeze(double(frames(x,y,:))),zRangeOrigin,'pchip');
        end
    end
    
%     h(1)=subplot(2,2,1);imagesc(squeeze(256*max(frames,[],1)/max(max(max(frames,[],1)))));
%     h(2)=subplot(2,2,2);imagesc(squeeze(256*max(newFrames,[],1)/max(max(max(newFrames,[],1)))));
%     linkaxes(h,'x','y');
    
    subplot(4,4,[1,2]);imagesc(zRange,yRange,squeeze(256*max(frames,[],1)/max(max(max(frames,[],1)))));%axis image;
    subplot(4,4,[5,6]);imagesc(zRangeOrigin,yRange,squeeze(256*max(newFrames,[],1)/max(max(max(newFrames,[],1)))));%axis image;
    
    subplot(4,4,[9,10]);imagesc(zRange,xRange,squeeze(256*max(frames,[],2)/max(max(max(frames,[],2)))));
    subplot(4,4,[13,14]);imagesc(zRangeOrigin,xRange,squeeze(256*max(newFrames,[],2)/max(max(max(newFrames,[],2)))));
        
    subplot(4,4,[3,4,7,8,11,12,15,16]);imagesc(yRange,xRange,squeeze(256*max(frames,[],3)/max(max(max(frames,[],3)))));axis image;   
end

delete(vid);
