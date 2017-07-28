
% fclose(esp_port);
% delete(esp_port);


% %% Setting up initial parameters for the com port
esp_port = serial('COM11');
esp_port.Baudrate=19200;
%esp_port.StopBits=1;
esp_port.Parity='none';
esp_port.Terminator='CR/LF';

%% Open port
fclose(esp_port);
fopen(esp_port);
status = esp_port.Status;

%% Sample power scan
sample = 'B2.1um';
wavelength = '960nm';
%% Set position to zero
% position_set = 0.0;
% smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(position_set)); fprintf(esp_port,smcErrCmnd); % set absolute position
% count = 0;
% smcErrCmnd = strcat(num2str(smcAxis),'TP'); smcErr = str2double(query(esp_port,smcErrCmnd));
% while (abs(smcErr-position_set)>0.01) && (count <100)
%     smcErr = str2double(query(esp_port,smcErrCmnd));
%     count =count+1;
%     pause(0.1);
% end
% fprintf('Rotation is set');

%% load calibration file
load('PowerCalibration960nm_20170712.mat');
figure;
plot(rotation,power,'o');
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);
hold on;
plot(f1)

title('Power calibration','FontSize', 20);
xlabel('Rotation (degrees)','FontSize', 20);
ylabel('Power (W)','FontSize', 20);
% set(gca,'FontSize',20);

p_max = f1.offset + f1.A;
p_min = f1.offset - f1.A;

%% Set rotation stage
smcAxis =1;
%
smcErrCmnd = strcat(num2str(smcAxis),'ID?');
smcErr = query(esp_port,smcErrCmnd);

%% Set to max power
% position_set = rot_calibr(p_max,f1);
% smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(position_set)); fprintf(esp_port,smcErrCmnd); % set absolute position
% count = 0;
% smcErrCmnd = strcat(num2str(smcAxis),'TP'); smcErr = str2double(query(esp_port,smcErrCmnd));
% while (abs(smcErr-position_set)>0.01) && (count <100)
%     smcErr = str2double(query(esp_port,smcErrCmnd));
%     count =count+1;
%     pause(0.1);
% end

%% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
src = getselectedsource(vid);
triggerconfig(vid,'manual');
ROIPosition = [768 768 512 512];
vid.ROIPosition = ROIPosition;
%src.ExposureTime = 0.01/2048*pixels; %in seconds
src.ExposureTime = 0.05;
vid.TriggerRepeat=0;
%     vid.FramesPerTrigger=1/scanningFrequency/src.ExposureTi me*nrCycles;
vid.FramesPerTrigger = 1;
folderName = 'E:\07-21 PowerScannin with B2.1um';


%% Scan power and acquire images

%p_var = p_min:0.05:p_max; %% Laser power measurement points
p_var = p_min:0.05:0.5; %% Laser power measurement points
p_var = p_min:0.05:0.5; %% Laser power measurement points

s1 = size(p_var);
clear im_array;
n_measurements = 10; % number of measurements per power point
%preview(vid);
figure(21);

vid.TriggerRepeat=s1(2);
start(vid);

for i=1:s1(2)
    
    position_set = rot_calibr(p_var(i),f1);
    smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(position_set)); fprintf(esp_port,smcErrCmnd); % set absolute position
    smcErrCmnd = strcat(num2str(smcAxis),'TP'); smcErr = str2double(query(esp_port,smcErrCmnd));
    count = 0;
    while (abs(smcErr-position_set)>0.01) && (count <1000)
        smcErr = str2double(query(esp_port,smcErrCmnd));
        count =count+1;
        pause(0.1);
    end
    pause(0.1);
    for j=1:n_measurements
        frame = getsnapshot(vid);
        im_array(:,:,(i-1)*n_measurements+j) = frame(:,:);
    end
    fprintf('Power: %d W\n',p_var(i));
    averageImage = mean(im_array(:,:,(i-1)*n_measurements+1:i*n_measurements),3);
    imagesc(averageImage);colorbar;
    fileName = strcat(folderName,'\image',num2str(i));
    %     imwrite(uint16(averageImage),strcat(fileName,'.tiff'),'tiff');
    
end
stop(vid);

position_set = rot_calibr(p_min,f1);
smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(position_set)); fprintf(esp_port,smcErrCmnd); % set

%% Save Data

% Working directory
today = datevec(date);
workdir = sprintf('E:/LightSheet/%d-%02d-%02d/',today(1),today(2),today(3));
if ~exist(workdir, 'dir')
    mkdir(workdir);
end
file_tiff = sprintf('%sPowerScan_%s_%s.tif', workdir,sample, wavelength);
% max_val = double(max(im_array(:)));
% im_array_norm = abs(double(im_array) ./ max_val);
npoints = size(im_array,3);
imwrite(im_array(:,:,1), file_tiff);
for i=1:npoints-1
    imwrite(im_array(:,:,i+1), file_tiff,'WriteMode','append');
end

% fprintf('block laser for taking background image');
% pause;
% for j=1:n_measurements
%     frame = getsnapshot(vid);
%     im_array(:,:,(i)*n_measurements+j) = frame(:,:);
% end

%stoppreview(vid);
%% Pre analyse data
cropx = [50 150];
%cropx = [250 300];
cropy = [150 300];



clear total_crop;
total_crop = squeeze(mean(mean(im_array(cropy(1):cropy(2),cropx(1):cropx(2),:),1),2));
profile = mean(reshape(total_crop, n_measurements, length(total_crop)/n_measurements),1);
% plot(log(p_var), log(profile),'ro');

% ft2 = fittype('a*x+b','independent',{'x'}, 'coefficients',{'a','b'});
% f2 = fit(log(p_var(10:end))',log(profile(10:end))',ft2);
% hold on;
% plot(f2)
%
% f2.a
%
%
% clear total_crop;
% total_crop = squeeze(mean(mean(im_array(cropy(1):cropy(2),cropx(1):cropx(2),:),1),2));
% figure;
% profile = mean(reshape(total_crop, n_measurements, length(total_crop)/n_measurements),1);
% plot(p_var, profile,'ro');
% ft3 = fittype('a*x^b+c','independent',{'x'}, 'coefficients',{'a','b','c'});
% f3 = fit(p_var',profile',ft3,'StartPoint',[1,2,0]);
% hold on;
% plot(f3)


profile_std = std(reshape(total_crop, n_measurements, length(total_crop)/n_measurements),1);
y_error = profile_std./ profile;

profile = profile - profile(1);
%profile =profile(2:end);
%y_error = y_error(2:end);
[m m_index]=max(profile(:));

figure;
subplot(1,2,1);
imagesc(mean(im_array(:,:,(m_index-1)*n_measurements+1:m_index*n_measurements),3));
hold on
rectangle('position',[cropx(1) cropy(1) cropx(2)-cropx(1) cropy(2)-cropy(1)])

title('Beam profile at max intensity','FontSize', 20);
xlabel('X (pixels)','FontSize', 20);
ylabel('Y (pixels)','FontSize', 20);
% Power calibration for different wavelengths;


% pmx = 1.09; pmn = 0.07; %% 900 nm
% pmx = 1.17; pmn = 0.06; %% 880 nm
% pmx = 1.35; pmn = 0.04; %% 860 nm
% pmx = 1.53; pmn = 0.03; %% 840 nm
% pmx = 1.52; pmn = 0.02; %% 820 nm
% pmx = 1.57; pmn = 0.01; %% 800 nm
% pmx = 1.63; pmn = 0.01; %% 780 nm
% pmx = 1.40; pmn = 0.01; %% 760 nm
% pmx = 1.07; pmn = 0.02; %% 740 nn
% p_var = p_var/pmx;
%
subplot(1,2,2);
% errorbar(log(p_var),real(log(profile(1:end-1))),real(y_error(1:end-1)),'*b');
% fitrange = [15 27]; %% range of experimental points included in the fit
% result1 = polyfit(log(p_var(fitrange(1):fitrange(2))),log(profile(fitrange(1):fitrange(2))),1);
% hold on; plot(log(p_var),polyval(result1,log(p_var)));
fitrange = [2 5]; %% range of experimental points included in the fit

result1 = polyfit(log(p_var(fitrange(1):fitrange(2))),log(profile(fitrange(1):fitrange(2))),1);
hold on; plot(log(p_var),polyval(result1,log(p_var)));

plot(log(p_var(fitrange(1):fitrange(2))),log(profile(fitrange(1):fitrange(2))),'ro');

formula = 'y = n*x + c';
coeff = sprintf('n = %.2f',result1(1));
text(0,3,coeff, 'FontSize', 20, 'FontWeight', 'bold','Color','blue');
text(-3,10,formula, 'FontSize', 20, 'FontWeight', 'bold','Color','black');

tit1 = sprintf('%s %s excitation', sample, wavelength);

title(tit1,'FontSize', 20);
xlabel('Log(Power[W])','FontSize', 20);
ylabel('Log(Fluorescence)','FontSize', 20);
set(gca,'FontSize',20);
hold off;
%plot(p_var,profile(1:end-1));


coeff
fclose(esp_port);

%% Plot all images in a single figure
% n_points = s1(2);
% fig = figure;
% col = floor(sqrt(n_points));
% row = floor(n_points/col) +1;
% scale = [0 max(im_array_norm(:))];
% %clims = [ 0 max_val];
% for i=1:n_points
%     fig = subplot(col,col+2,i);
%     imagesc(mean(im_array(:,:,(i-1)*n_measurements+1:i*n_measurements),3),scale);
%     set(gca,'visible','off')
%     axis equal;
% end
%colormap(gray);