%% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
src = getselectedsource(vid);
triggerconfig(vid,'manual');
pixelNumber = 1024;
vid.ROIPosition = [(2048-pixelNumber)/2 (2048-pixelNumber)/2 pixelNumber pixelNumber];
%src.ExposureTime = 0.01/2048*pixels; %in seconds
src.ExposureTime = 1;
vid.TriggerRepeat=0;
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTime*nrCycles;
vid.FramesPerTrigger = 1;

% each trigger will include 2 full cycle, i.e. four volumes of images


FOV = 580;
pixelSize = 580/2048;
h = figure;

%         start(vid);

for i=1:100
    
    img = getsnapshot(vid);
    % xlsread(coordinateFileName);
    imgCenter = uint16(zeros(pixelNumber,pixelNumber));
    imgCenter(50:(pixelNumber-50),50:(pixelNumber-50))=img(50:(pixelNumber-50),50:(pixelNumber-50));
    beadCoordinate = FastPeakFind((imgCenter),0.2*max(imgCenter(:)));
    resultNumber = size(beadCoordinate);
    beadCoordinate = reshape(beadCoordinate, 2, resultNumber(1)/2);
    
    for i = 1:length(beadCoordinate(2,:))
        %get the image of the bead based on the coordinate 11*11 pixels
        bead = double(img((beadCoordinate(2,i)-5):(beadCoordinate(2,i)+5),(beadCoordinate(1,i)-5):(beadCoordinate(1,i)+5)));
        % plot the profile if the beads
        profile = sum(bead,1);
        % interpolate 5 times, result in 51*51 pixels;
        profile2 = interp1(1:11,profile,1:0.2:11,'pchip');
        % normalize to 0~1
        profile2 = (profile2-min(profile2))/(max(profile2)-min(profile2));
        % find the number of pixels above 0.5, as the indicator of the
        % focus quality, (the smaller, the better)
        beadCoordinate(3,i) = sum(profile2>0.5);
        % save the profile of the bead to beadProfiles
        beadProfiles(1:51,i) = profile2;
        %        subplot(2,1,1);plot(profile);xlim([0 11]);
        %         subplot(2,1,2);plot(profile2);
    end
    
    %sort the beads according to the quality of focus (numbers above 0.5)
    [sorted_x,index] = sort(beadCoordinate(3,:),'ascend');
    %average numbers of beads profile;
    beadNumber = 20;
    sumProfile = sum(beadProfiles(:,index(1:beadNumber)),2)/beadNumber;
    
    % fit the sumProfile to a Gaussian function
    x=pixelSize/5:pixelSize/5:51*pixelSize/5;
    x = x-51*pixelSize/5/2;
    f = fit(x',sumProfile,'gauss1');
    resolution = 2*f.c1*sqrt(log(2));
    
    %     subplot(1,2,1),
    imagesc(img);axis image;
    hold on;
    plot(beadCoordinate(1,:),beadCoordinate(2,:),'r*');axis image;
    hold off;
%     subplot(1,2,2),plot(f,x,sumProfile);
    %     xlabel('\mum');ylabel('Normalized intensity');
        str = sprintf('FWHM is %.2f ',resolution);
        title(strcat(str,'\mum'));
    
    %     set(findall(gcf,'type','text'),'FontSize',25)
    %     AX=legend('Raw data','Fitted curve');
    %     AX.FontSize = 20;
    
end
stop(vid);

delete(vid);

%     resolution2 = 2.355*f.c1/sqrt(2);