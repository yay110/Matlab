%% One has to initialize and set up Andor camera and newportcontroller before using this script
%% Choose initial stage settings
depth_low = -0.05; % in mm
depth_high = 0.05;
depth_step = 0.001;
depth_set = -0.1;

%% initialize stage position
smcAxis =1; %% controller axis
smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(depth_set)); fprintf(esp_port,smcErrCmnd);
smcErrCmnd = '1TP'; smcErr = str2double(query(esp_port,smcErrCmnd)) %% check absolute position

count = 0;
while (abs(smcErr-depth_set)>0.0001) && (count <10)
    smcErr = str2double(query(esp_port,smcErrCmnd));
    count =count+1;
    pause(0.001);
end

%% set up camera settings

%% Get detector size
i_x = 0; i_x_pt = libpointer('int32Ptr', i_x);
i_y = 0; i_y_pt = libpointer('int32Ptr', i_y);
[err_code, i_x, i_y] = calllib('atmcd64d','GetDetector', i_x_pt, i_y_pt);  %% ui_error = GetDetector(&i_x,&i_y);
%% Binning
im_binning = 4; max_x = 4*floor(double(i_x)/im_binning); max_y = 4*floor(double(i_y)/im_binning);
err_code = calllib('atmcd64d','SetImage',im_binning,im_binning,1,max_x,1,max_y);
readouttime = 0; readouttime_pt = libpointer('singlePtr', readouttime);
[err_code, readouttime_pt] = calllib('atmcd64d','GetReadOutTime',readouttime_pt);

i_x = floor(double(i_x)/im_binning); i_y = floor(double(i_y)/im_binning);
im_array = zeros(i_x, i_y); im_array_pt = libpointer('longPtr', im_array);
size_xy = i_x * i_y;
%% Exoposure
err_code = calllib('atmcd64d','SetExposureTime',0.1); %%   ui_error = SetExposureTime(params.f_exposureTime);


% tic
% err_code = calllib('atmcd64d','StartAcquisition');
% toc
% err_code = calllib('atmcd64d','WaitForAcquisition');
% toc
% [err_code, im_array] = calllib('atmcd64d','GetMostRecentImage',im_array_pt, size_xy);   %% ui_error = GetMostRecentImage(lp_data,ul_size);
% figure;
% imagesc(im_array, [0.0 5500]);  

%% Prepare profile analysis arrays
npoints = 100;
clear im_array_scan;
im_array_scan = zeros(i_x,i_y,npoints+1);
profile = zeros(1,npoints+1);
depth  = zeros(1,npoints+1);
x_min = floor(300/im_binning); x_max = floor(700/im_binning);
y_min = floor(300/im_binning); y_max = floor(700/im_binning);
step = 1; %% depth step in um
offset = 0;
%% Run septh scan and acqusition
for i=0:npoints
    
    depth_set = depth_low + depth_step*i;
    smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(depth_set)); fprintf(esp_port,smcErrCmnd);
    count = 0;
    smcErrCmnd = '1TP'; smcErr = str2double(query(esp_port,smcErrCmnd));
    while (abs(smcErr-depth_set)>0.0001) && (count <100)
        smcErr = str2double(query(esp_port,smcErrCmnd));
        count =count+1;
        pause(0.001);
    end
    pause(0.05);
    %tic
    err_code = calllib('atmcd64d','StartAcquisition');
    err_code = calllib('atmcd64d','WaitForAcquisition');
    [err_code, im_array] = calllib('atmcd64d','GetMostRecentImage',im_array_pt, size_xy);
    %toc
    im_array_scan(:,:,i+1) = im_array(:,:);
    %toc
    fprintf('Depth: %f \n',depth_set);
end

fprintf('block laser for taking background image');
pause;
err_code = calllib('atmcd64d','StartAcquisition');
err_code = calllib('atmcd64d','WaitForAcquisition');
[err_code, im_array] = calllib('atmcd64d','GetMostRecentImage',im_array_pt, size_xy);
figure; imagesc(im_array);
im_bkg = im_array;
%err_code = calllib('atmcd64d','CancelWait');
%bsxfun

%% Export images into tiff stack
% Working directory 
today = datevec(date);
workdir = sprintf('F:/TFproject/DepthScans/%d%02d%02d/',today(1),today(2),today(3));
if ~exist(workdir, 'dir')
  mkdir(workdir);
end
file_tiff = sprintf('%sDepthScan_step_%dnm.tif', workdir,int32((depth_step)*1000000));
max_val = max(max(max(im_array_scan)));
im_array_scan = abs(im_array_scan ./ max_val);
bkg_max = double(max(max(im_array)));
if max_val < bkg_max
    max_val = double(max(max(im_array)));
end
im_bkg = double(im_array) ./ max_val;



imwrite(im_array_scan(:,:,1), file_tiff);
for i=1:npoints  
    imwrite(im_array_scan(:,:,i+1), file_tiff,'WriteMode','append');
end
imwrite(im_bkg, file_tiff,'WriteMode','append');


%% Export images into video file 
% max_val = max(max(max(im_array_scan)));
% im_array_scan = abs(im_array_scan ./ max_val);
% bkg_max = double(max(max(im_array)));
% if max_val < bkg_max
%     max_val = double(max(max(im_array)));
% end
% im_bkg = double(im_array) ./ max_val;
% workdir = 'F:\TFproject\DepthScans\';
% file_avi = sprintf('%snewfile.avi', workdir);
% v = VideoWriter(file_avi,'Grayscale AVI');
% v.FrameRate = 1;
% open(v);
% writeVideo(v,im_array_scan);
% writeVideo(v,im_bkg);
% close(v);

%% Subtract the background and calculate the depth profile
%im_array_scan = im_array_scan - double(repmat(im_array,[1 1 npoints+1]));
for i=0:npoints    
    %im_array_scan(:,:,i+1) = im_array_scan(:,:,i+1) - im_array(:,:);
    profile(i+1) = sum(sum(im_array_scan(y_min:y_max,x_min:x_max,i+1)));    
    depth(i+1) = i*step + offset;
end
profile = profile - profile(end); %% if the last frame is the background
%% Analyse depth profile
norm_prof = profile/max(abs(profile(:)));
[m m_index]=max(norm_prof(:));

figure(3);
subplot(1,2,1);
hold on;
plot(depth-depth(m_index),norm_prof);

subplot(1,2,2); 
imagesc(im_array_scan(:,:,m_index));
hold on
rectangle('position',[y_min x_min y_max-y_min x_max-x_min])


smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(depth_set)); fprintf(esp_port,smcErrCmnd);
pause(2);
smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(0)); fprintf(esp_port,smcErrCmnd);
