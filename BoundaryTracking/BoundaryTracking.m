close all
clc
clear


%% Determining filenames of data to be analysed
cd('C:\Users\yay-1_000\OneDrive\General in St Andrews\projects\2015 Oxford- neuron cell paper\revision and related\');                %videos saved in portable harddrive
dataFolder        	= 'stability check';                         	%Subfolder
% FOLDER           	= [cd, filesep, dataFolder];           	%Full path to folder containing the data files
FOLDER              = dataFolder;
TYPE               	= 'avi';                               	%File format
SIZE              	= [10, 3000]*10^6;                    	%Range of files size in B
[FN]              	= filenames(FOLDER,TYPE,SIZE,false);    %List of filenames in subfolder with the format 'TYPE' and file size 'SIZE'.


%% Constants
p2m               	= 0.147*1e-6;                      %Pixel to meter [m/pix]
T               	= 297.15;                           %Temperature [K] 24 degrees
kB               	= 1.3806488e-23;                    %Boltzmann constant [m^2*kg/(s^2*K)]


%% Displaying parameters
portionToProcess    = 1;
ShowPreview         = 0;        %If display preview of tracking

for ID = 1:1%size(FN,1)
    tic
    
    %% Loading movie
    inputFile           = [FOLDER, filesep, FN{ID},'.', TYPE];
    inputVideo          = VideoReader(inputFile);
    
    % Variables
    nFrames             = portionToProcess * inputVideo.NumberOfFrames;        %Number of frames [#]
    fps                 = inputVideo.FrameRate;             %Camera framerate [Hz]
    fps                 = 800;             %Camera framerate [Hz]

    t                	= 0:1/fps:(nFrames-1)/fps;          %Time vector [s]
    
    % Interpolation factor
    F                	= 5;                                %Artificially increase resolution by a factor of F. F = 1 => no interpolation.
    
    %     % Determining normalisation constants
    %     I0mx(1:nFrames)=0;
    %     I0mn(1:nFrames)=0;
    %     for k = 1:nFrames
    %         I0                      = read(inputVideo, k);
    %         I0mx(k)                 = max(double(I0(:)));
    %         I0mn(k)                 = min(double(I0(:)));
    %
    %         if k/round(nFrames/100) == round(k/round(nFrames/100))
    %             clc
    %             display(['Determining normalisation constants...',sprintf('\nMovie ID: %1.0f', ID), sprintf('\nProgress: %1.1f%%', 100*k/nFrames)])
    %         end
    %     end
    %     I0max                   = max(I0mx);
    %     I0min                   = min(I0mn);
    I0max=255;
    I0min=1;
    I0                      = read(inputVideo, 1);
    [Xin,Yin]             	= meshgrid(1:size(I0,1));
    [Xout,Yout]             = meshgrid(1:1/F:size(I0,1));
    
    %% Boundary Tracking
    Br=cell(1,nFrames);
    Bc=cell(1,nFrames);
    r0(1:nFrames)=0;
    c0(1:nFrames)=0;
    for n = 1:nFrames
        %load frame
        I0                  = read(inputVideo, n);
        I0Norm              = (double(I0)-I0min)/(I0max-I0min);
        I                   = interp2(Xin,Yin,I0Norm(:,:,1),Xout,Yout,'linear');
        
        % Color to black and white
        mu              	= mean(I(:));
        sd               	= std(I(:));
%        th               	= 1.9*(mean(I(:))+min(I(:)))/2;             %Threshold
        th                  = 0.1;
        bw               	= im2bw(I,th);                          %Apply threshold to binary
        
        bw0                 = bw;%imcomplement(bw);                     % imcomplement the image
        % Cleaning up image
        n1               	= 10;
        bw1               	= bwareaopen(bw0,n1);                   %n1 = 100, get rid of any structure smaller than 100 pixels
        n2                  = 1;
        bw2              	= imclose(bw1,strel('disk', n2));       %n2 = 1, close the structure
        n3                  = 0;
        bw3                 = bwareaopen(bw2,n3);                   %n3 = 0, function not used
        bw4              	= imfill(bw3,'holes');
        
        % Get the boundary
        B                 	= bwboundaries(bw4,'noholes');
        Br{n}             	= B{1}(:,1);
        Bc{n}             	= B{1}(:,2);
        
        % Get the position
        r0(n)             	= mean(Br{n});
        c0(n)              	= mean(Bc{n});
        
        % Plotting
        if ShowPreview
            subplot(231), imagesc(I),   axis image, title('original image')
            subplot(232), imagesc(bw), axis image, title('binary image')
            subplot(233), imagesc(bw0), axis image, title('imcomplemented image')
            subplot(234), imagesc(bw1), axis image, title('cleared image')
            subplot(235), imagesc(bw2), axis image, title('closed image')
            subplot(236), imagesc(I), axis image, title('displaying center of image'), colormap(gray)
            hold on
            plot(Bc{n},Br{n},'r')
            plot(c0(n),r0(n),'rx','MarkerSize',10)
            drawnow
        end
        
        % Progress
        if n/round(nFrames/100) == round(n/round(nFrames/100))
            clc
            display(['Tracking...',sprintf('\nMovie ID: %1.0f', ID), sprintf('\nProgress: %1.1f%%', 100*n/nFrames)])
        end
    end
    
    % Centre data at zero and convert unit to meter.
    r(:,ID)             = p2m*(r0 - mean(r0))./F;
    c(:,ID)         	= p2m*(c0 - mean(c0))./F;
    
    % Get the standard deviation
    rstd(ID)            = std(r(:,ID));
    cstd(ID)            = std(c(:,ID));
    
    % Get the force constant
    rkappa(ID)         	= 1e6*kB*T./var(r(:,ID));                     %[pN/um]
    ckappa(ID)          = 1e6*kB*T./var(c(:,ID));                     %[pN/um]
    
    Time(ID)            = toc;                                  %Computing time per movie
    
    
    %% Plotting
    scrnsz              = get(0,'ScreenSize');
    m2um                = 1e6;                                  %meter to micrometer
    wMax                = m2um*max([max(abs(r(:,ID))),max(abs(c(:,ID)))]);
    wLim                = 1.1*[-wMax, wMax];
    
    nhist               = 50;
    h = figure('Position', scrnsz);
    subplot(2,2,1),plot(t,m2um*r(:,ID))
    title('X direction Positional Trace'),xlabel('time [s]'),ylabel('position [um]'),xlim([t(1),t(end)]),ylim(wLim)
    subplot(2,2,2),plot(t,m2um*c(:,ID))
    title('Y direction Positional Trace'),xlabel('time [s]'),ylabel('position [um]'),xlim([t(1),t(end)]),ylim(wLim)
    subplot(2,2,3),hist(m2um*r(:,ID),nhist)
    title(sprintf('Histogram, std = %1.2f um',m2um*rstd(ID))),xlabel('position [um]'),ylabel('Counts [#]'),xlim(wLim)
    subplot(2,2,4),hist(m2um*c(:,ID),nhist)
    title(sprintf('Histogram, std = %1.2f um',m2um*cstd(ID))),xlabel('position [um]'),ylabel('Counts [#]'),xlim(wLim)
     saveas(h,FN{ID},'fig')
     saveas(h,FN{ID},'jpg')
%    close all;
    X=r(:,ID);
    Z=c(:,ID);
    save  ((strcat('stiffnessFeb04-',FN{ID},'.mat')), 'X', 'Z');
end
%save (strcat('stiffnessFeb04',FN{ID},'.mat'))

save stiffnessFeb04-1.mat
