%Give the folder name and the file name of the images, it then go through
%the beads, find the best focucsed ones, average them and do Gaussian
%fitting.

% Created by Zhengyi Yang on 29th August 2016.

folderName = 'f:\Acoustic trapping paper related\resolution optimization';

for n = [46,37,38,40]%[36:43,46,47]
    FOV = 580;
    pixelSize = 580/2048;
    
    imageFileName = strcat(folderName, '\Image000',num2str(n),'.tif');
    % coordinateFileName = strcat(folderName,'\Results',num2str(n),'.xls');
    img = imread(imageFileName);
    % xlsread(coordinateFileName);
    imgCenter = uint16(zeros(2048,2048));
    imgCenter(100:2000,100:2000)=img(100:2000,100:2000);
    beadCoordinate = FastPeakFind((imgCenter),0.05*max(imgCenter(:)));
    resultNumber = size(beadCoordinate);
    beadCoordinate = reshape(beadCoordinate, 2, resultNumber(1)/2);
    
    %     figure;
    %     h(1) = subplot(1,2,1);imagesc(img);axis image;
    %     h(2) = subplot(1,2,2);imagesc(img);axis image;
    %     hold;
    %     plot(beadCoordinate(1,:),beadCoordinate(2,:),'r*');axis image;
    %     plot(108,580,'ro')
    %
    %     hold;
    %     linkaxes(h);
    
    
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
    resolution(n) = 2*f.c1*sqrt(log(2));
    
    h = figure;
    plot(f,x,sumProfile);
    xlabel('\mum');ylabel('Normalized intensity');
    str = sprintf('FWHM is %.2f ',resolution(n));
    title(strcat(str,'\mum'));
    %     title(['The resolution is ', num2str(resolution(n)),' \mum']);
    
    %      set(0,'DefaultAxesFontSize', 14)
    %      axes('FontSize',24);
    set(findall(gcf,'type','text'),'FontSize',25)
    figureName = strcat(folderName, '\Image000',num2str(n));
    AX=legend('Raw data','Fitted curve');
    AX.FontSize = 20;
    saveas(h,figureName,'svg');
    saveas(h,figureName,'png');
    
    
    %     resolution2 = 2.355*f.c1/sqrt(2);
end