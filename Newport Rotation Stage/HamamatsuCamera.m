%vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode');
vid = videoinput('hamamatsu', 1, 'MONO16_BIN2x2_1024x1024_FastMode');

src = getselectedsource(vid);
%get(src); % Get properties of the camera

%vid.FramesPerTrigger = 1;
%%
exposure = 0.6;
src.ExposureTime = exposure;
%src.SensorCooler = 'ON';


%%
%figure(11);
%colormap(gray)
preview(vid);
%%
% triggerconfig(vid, 'manual');
%set(vid, 'TriggerRepeat', Inf);
% start(vid);
npoints = 10000;
figure(20);
for i=0:npoints
    %pause(exposure);
    frame = getsnapshot(vid);
    %curFrame = getdata(vid);
    imagesc(frame);
    %imagesc(frame,[0 1000] );
    %colormap(gray);
    colorbar;
    axis image;
    drawnow;
    kkey = get(gcf,'CurrentCharacter');
    if (kkey)
        break;
    end
end
% stop(vid);
%%
stoppreview(vid);


delete(vid);