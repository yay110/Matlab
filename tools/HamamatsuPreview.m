%% camera initialisation
% vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
vid = videoinput('hamamatsu', 1, 'MONO16_BIN2x2_1024x1024_FastMode');
%  vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode');
src = getselectedsource(vid);
triggerconfig(vid,'manual');

% ROIPosition = [0 0 2048 2048];

% vid.ROIPosition = ROIPosition;
%src.ExposureTime = 0.01/2048*pixels; %in seconds
src.ExposureTime = 0.1;
vid.TriggerRepeat=10000;
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
vid.FramesPerTrigger = 1;
folderName = 'E:\01-31 PowerScannin with fluorescein';
start(vid);

gca = figure(1);
axis image;
set(gca,'Units','pixels')
set(gca,'Position',[100 100 900 800]);
pos = get(gca, 'Position');
trigger(vid);
frame = getdata(vid);
imagesc(frame);
axis image;
colorbar;
colormap(gray);
set(gcf,'doublebuffer','off');

%     set(gca, 'xlimmode','manual',...
%            'ylimmode','manual',...
%            'zlimmode','manual',...
%            'climmode','manual',...
%            'alimmode','manual');

for i = 0:20000
    %     tic
    trigger(vid);
    
    %     frame = getsnapshot(vid);
    frame = getdata(vid);
%     frame = insertMarker(frame,[512 512]);
    imagesc(frame);
    axis image;
    maxV = max(frame(:));
    minV = min(frame(:));
    title(strcat('Dynamic Range is from Min=',num2str(floor(minV/4)),' to Max=',num2str(ceil(maxV/4))));
    %     colorbar;
    %     set(gcf,'currentchar',' ')         % set a dummy character
    %     while get(gcf,'currentchar')=='c'  % which gets changed when key is pressed
    %
    %         break
    %
    %     end
    %     toc
end

stop(vid);
