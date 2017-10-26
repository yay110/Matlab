clear;
close all;
tic;
% %% Setting up initial parameters for the com port
stage = NewportESP();
laser = Chameleon();
cd 'C:\Users\OMG\Documents\GitHub\Matlab\Newport Rotation Stage';
%% camera initialisation
Andor = Andor();
[ret,temperature] = GetTemperature();
fprintf('Sensor Temperature is %d\n',temperature);
Andor.AcquisitionMode = 1;
Andor.ExposureTime = 0.2;
Andor.NumberKinetics = 10;
cropSize = 40;


%% Sample power scan
sample = 'testing';

% Working directory
today = datevec(date);
folderName = sprintf('D:/LightSheet/%d-%02d-%02d PowerScann in BBO 1E-4M/',today(1),today(2),today(3));
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

p_step = 0.015;                  %% Laser power measurement points
p_scan = 0.3;
wavelengthRange = 800;%:50:1050;
fig = figure('units','normalized','outerposition',[0 0 1 1]);

for ii = 1:length(wavelengthRange)
    wavelength = wavelengthRange(ii);
    laser.Wavelength = wavelength;
    %% load calibration file and fit into sin wave
    file = dir(['*doubleWP-PowerCalibration' num2str(wavelength) 'nm*.mat']);
    load(file.name);
    
    ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
    f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);
    p_max = f1.offset + f1.A;
    p_min = f1.offset - f1.A;
    
    if p_max < p_scan
        p_scan = p_max;
        p_step = p_max/10;
    end
    p_var = p_min:p_step:p_scan;                  %% Laser power measurement points
    position_set = rot_calibr(p_min,f1);
    stage.absPosition = position_set;
    smcErr = stage.absPosition;
    count = 0;
    while (abs(smcErr-position_set)>0.01) && (count <1000)
        smcErr = stage.absPosition;
        count =count+1;
        pause(0.1);
    end
    pause(0.1);
    
    %% Scan power and acquire images
    fitrange = 2:length(p_var); %% range of experimental points included in the fit
    
    s1 = size(p_var);
    clear im_array;
    n_measurements = 10; % number of measurements per power point
    subplot(2,2,1)
    
    
    %         for j=1:n_measurements
    %             frame = Andor.getImage;
    %             background(:,:,j) = int16(frame);
    %         end
    %         fprintf('taking background image;');
    %         averageBackground = mean(background(:,:,:),3);
    %                 averageBackground(:,1:50) = 500;
    %         averageBackground(end-50:end,:) = 500;
    %         imagesc(averageBackground);
    %         colorbar;
    %         tit1 = sprintf('background');
    %
    %         title(tit1,'FontSize',10);        fileName = strcat(folderName,'\image',num2str(wavelength));
    %         label = sprintf('maximum in the image is %s', num2str(max(averageBackground(:))));
    %         xlabel(label);
    
    
    for i=1:s1(2)
        if i==2
            laser.SOpen;
        end
        position_set = rot_calibr(p_var(i),f1);
        stage.absPosition = position_set;
        smcErr = stage.absPosition;
        count = 0;
        while (abs(smcErr-position_set)>0.01) && (count <1000)
            smcErr = stage.absPosition;
            count =count+1;
            pause(0.1);
        end
        pause(0.1);
        for j=1:n_measurements
            frame = Andor.getImage;
            im_array(:,:,(i-1)*n_measurements+j) = int16(frame);
        end
        fprintf('Power: %d W\n',p_var(i));
        averageImage = mean(im_array(:,:,(i-1)*n_measurements+1:i*n_measurements),3);
        averageImage(:,1:50) = 500;
        averageImage(end-50:end,:) = 500;
        imagesc(averageImage);
        colorbar;
        tit1 = sprintf('%s %snm excitation', sample, num2str(wavelength));
        
        title(tit1,'FontSize',10);        fileName = strcat(folderName,'\image',num2str(wavelength));
        label = sprintf('maximum in the image is %s', num2str(max(averageImage(:))));
        xlabel(label);
        %         imwrite(uint16(averageImage),strcat(fileName,'.tiff'),'tiff');
    end
    
    laser.SClose;
    
    
    % Save Data
    
    % Working directory
    today = datevec(date);
    file_tiff = sprintf('%sPowerScan_%s_%s.tiff', folderName,sample, num2str(wavelength));
    % max_val = double(max(im_array(:)));
    % im_array_norm = abs(double(im_array) ./ max_val);
    npoints = size(im_array,3);
    imwrite(uint16(im_array(:,:,1)), file_tiff);
    for i=1:npoints-1
%         pause(0.1);
        imwrite(uint16(im_array(:,:,i+1)), file_tiff,'WriteMode','append');
    end
    %     averageImage = medfilt2(averageImage);  % remove salt and pepper noise
    if  max(averageImage(:))< 500
        fprintf('There is no signal!\n')
    else
        
        %         if ii==1
        [centerX,centerY]= find(averageImage == max(averageImage(:)));
        %stoppreview(vid);
        %% Pre analyse data
        cropx = [centerX(1)-cropSize/2 centerX(1)+cropSize/2];
        %cropx = [250 300];
        cropy = [centerY(1)-cropSize/2 centerY(1)+cropSize/2];
        %         end;
        
        subplot(2,2,3)
        imagesc(averageImage(cropx(1):cropx(2),cropy(1):cropy(2)));
        clear total_crop;
        total_crop = squeeze(mean(mean(im_array(cropx(1):cropx(2),cropy(1):cropy(2),:),1),2));
        profile = mean(reshape(total_crop, n_measurements, length(total_crop)/n_measurements),1);
        
        
        profile_std = std(reshape(total_crop, n_measurements, length(total_crop)/n_measurements),1);
        y_error = profile_std./ profile;
        
        profile = profile - profile(1);
        [m,m_index]=max(profile(:));
        
        % figure;
        % subplot(1,2,1);
        % imagesc(mean(im_array(:,:,(m_index-1)*n_measurements+1:m_index*n_measurements),3));
        hold on
        %     rectangle('position',[cropx(1) cropy(1) cropx(2)-cropx(1) cropy(2)-cropy(1)])
        %
        % title('Beam profile at max intensity','FontSize', 20);
        % xlabel('X (pixels)','FontSize', 20);
        % ylabel('Y (pixels)','FontSize', 20);
        % Power calibration for different wavelengths;
        subplot(2,2,2)
        plot(p_var,profile);
        subplot(2,2,4);
        
        result1 = polyfit(log(p_var(fitrange)),log(profile(fitrange)),1);
        plot(log(p_var(2:end)),polyval(result1,log(p_var(2:end))));
        hold on;
        % plot(log(p_var),log(profile));
        
        plot(log(p_var(fitrange)),log(profile(fitrange)),'ro');
        hold off;
        
        formula = 'y = n*x + c';
        coeff = sprintf('n = %.2f',result1(1));
        text(0,3,coeff, 'FontSize', 10, 'FontWeight', 'bold','Color','blue');
        text(-3,10,formula, 'FontSize', 10, 'FontWeight', 'bold','Color','black');
        
        tit1 = sprintf('%s %snm excitation n = %.2f', sample, num2str(wavelength),result1(1));
        
        title(tit1,'FontSize',10);
        xlabel('Log(Power[W])','FontSize', 10);
        ylabel('Log(Fluorescence)','FontSize', 10);
        set(gca,'FontSize',10);
        hold off;
        % plot(p_var,profile(1:end-1));
        saveas(fig,strcat(folderName,'\',tit1,'.png'))
        saveas(fig,strcat(folderName,'\',tit1,'.fig'))
        
        coefficient(ii) = result1(1);
    end
end
laser.Wavelength = 950;
% stage.absPosition =  rot_calibr(p_min,f1);
% stage.absPosition =  rot_calibr(0.2,f1);

laser.Close;
stage.Close;
coefficient
runTime = toc/60;
fprintf('Test took %s mins\n',runTime);