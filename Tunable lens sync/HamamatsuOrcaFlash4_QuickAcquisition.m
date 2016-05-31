%%% Initialise camera
clear frames;
disp('Initialising the Hamamatsu camera...')
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
%vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode');

triggerconfig(vid,'manual')
src = getselectedsource(vid);
vid.FramesPerTrigger=1;
vid.TriggerRepeat=Inf;
%vid.ROIPosition = [0 0 512 512];

start(vid)
pause(1)

%%%

src.ExposureTime = 0.01; %in seconds
frameNumber = 50;

pixels = 2048;
pixelH = pixels;
pixelV = pixels;

%% Preview camera
% trigger(vid);
% previewFrame=getdata(vid,1);
% previewFig=figure();
% imagesc(previewFrame);axis image;
% %[~,~]=ginput(1);
% close(previewFig)

%%%

frames = uint16(zeros(pixelH,pixelV,frameNumber));
times=zeros(1,frameNumber);

disp('Hamamatsu camera initialisation complete')


for n =1:frameNumber
    tic
    trigger(vid);
    % [frames(:,:,n),times(n),absTimes(n)] = getdata(vid);

    pause(0.01);%     frames(:,:,n) = getsnapshot(vid);
        times(n) = toc;
end
tic;
frames = getdata(vid);
toc;
%  for n=1:frameNumber
%      tic
%     frames(:,:,n)=getdata(vid,1);
% %     times(n)=toc;
%     % Acquire frame
%      times(n)=toc;
% 
% 
%     %%% Live camera display
%     % figure(liveDisplayFig);
%     %        if mod(n,10) == 0
%     % %        imagesc(frames(:,:,n));%axis image;
%     %        end
%     %       drawnow;shg;
%     %%%
%  end

figure(55);
% imagesc(frames(:,:,1));
plot(times);
delete(vid);