%%% Initialise camera
    disp('Initialising camera...')
    vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode');
    triggerconfig(vid,'manual')
    src = getselectedsource(vid);
    src.ExposureTime = 0.02; %in seconds (I think)
    vid.FramesPerTrigger=1;
    vid.TriggerRepeat=Inf;
    start(vid)
    pause(1)
    disp('Camera initialisation complete')
    %%%
    
    disp('Setting other parameters...')
    %%% Preview camera
    % Acquire frame
    previewFrame=zeros(frameV,frameH);
    trigger(vid);
    previewFrame=getdata(vid,1);
    previewFig=figure();imagesc(previewFrame);axis image;
    title('Click on point in image to perform correction on to continue')
    [~,~]=ginput(1);
    close(previewFig)
    %%%
    
    for n=1:1000
        % Acquire frame
        trigger(vid);
        frames(:,:,n)=getdata(vid,1);
        %%% Live camera display
        figure(liveDisplayFig);
        imagesc(frames(:,:,n));axis image;
        drawnow;shg;
        %%%
    end