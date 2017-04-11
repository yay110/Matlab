clear;

folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\2nd sample from 1\0. control without drug';
folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\1. before drug';
folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\4. after drug';
folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\5. after drug. take 2';
%   folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\9. Trippling drug';
%   folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\10. after trippling drug';
folderName = 'F:\2017-03-08 zebrafish heart beat\zebrafish Mar08\3rd sample from 0\11. 3nd sample long term after trippling the drug';
folderName = 'E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\2.adding and after drug';

folderNames = dir(folderName);

if folderNames(3).isdir
    for n = 1:length(folderNames)-2
        subFolderNames(n).name = strcat(folderName,'\',folderNames(n+2).name);
    end
    sectionLength = 9;
else 
    subFolderNames.name = folderName;
    sectionLength = 9;
end

for nn = 34%1:length(subFolderNames)
    
    cd(subFolderNames(1).name)
    fileName = dir('**.tiff');
    
    Fs = 20;            %sampling frequency
    frameNumber = Fs*sectionLength;
    %read the images
    
    sectionNumber = floor(length(fileName)/(frameNumber));
    
    for n=1:sectionNumber
        for m = 1:frameNumber
            subDataCube(:,:,m) = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+m).name));
        end
        subDataCube = subDataCube-100;
        subDataCube(subDataCube<0) = 0;
        subDataCube = double(subDataCube)./double(max(subDataCube(:)));
        
        for m = 1:frameNumber
            I = subDataCube(:,:,m);
            
            BW = im2bw(I, 0.005);
                         subplot(1,2,1),imagesc(I);
                         subplot(1,2,2),imagesc(BW);
            n1               	= 200;
            BW1               	= bwareaopen(BW,n1);
            n2                  = 1;
            BW2              	= imclose(BW1,strel('disk', n2));       %n2 = 1, close the structure
            BW3              	= imfill(BW2,'holes');
                         imagesc(BW3)
                         hold on;
            
            B = bwboundaries(BW3,'noholes');
                         for k = 1:length(B)
                           boundary = B{k};
                            plot(boundary(:,2), boundary(:,1));
                        end
                        hold off;
            
            heartBoundary(n,m) = bwarea(BW3);
            %             measurements = regionprops(BW3, 'Area');
            %             heatBoundary(m) = [measurements.Area];
        end
        n
    end
    
    
    %     cubeIdx = 1;
    %     for n = 1:length(fileNames)
    %         if (~fileNames(n).isdir) %&& contains(fileNames(n).name,format)
    %             %             dataCube(:,:,cubeIdx) = imread(subFolderNames(nn).name,'\',fileNames(n).name));
    %             frame = imread(strcat(subFolderNames(nn).name,'\',fileNames(n).name));
    %             profile(cubeIdx) = sum(sum(frame()));
    %             cubeIdx = cubeIdx +1;
    %         end
    %     end
    %
    
    %     dataCube = readImagesFromFolder(subFolderNames(nn).name);
    %
    %     %find maxima spot
    %     firstFrame = dataCube(:,:,1);
    %     [row, col] = find(firstFrame == max(firstFrame(:)));
    %
    %     % row in matlab is Y in imagej
    %     % col in matlab is X in imagej
    %
    %     %Z profile near area of maxima
    %     %         dataCube = dataCube(:,:,1:1800);
    %
    %     for n = 1:size(dataCube,3)
    %         %             profile(n) = sum(sum(dataCube(row-20:row+20, col-20:col+20,n)));
    %         profile(n) = sum(sum(dataCube(:,:,n)));
    %
    %     end
    
    
    %define parameters of the signal for FFT
    Fs = 20;                        %sampling frequency
    T = 1/Fs;                       % sampling period
    L = sectionLength*Fs;       % Length of signal
    t = (0:L-1)*T;                  %Time vector
    
    %resize the profile to small sections
    %          sectionLength = 10;         % Unit second to divide the images
    %     sectionFrames = sectionLength*Fs;
    
    
    %     profile = profile(1:floor(length(profile)/(sectionLength*Fs))*(sectionLength*Fs));
    %     sectionProfile = reshape(profile, [sectionFrames length(profile)/sectionFrames]);
    
    % plot(t,profile)
    % % title('Signal Corrupted with Zero-Mean Random Noise')
    % xlabel('t (seconds)')
    % ylabel('Profile')
    clear beatFreqency;
    % figure;
    % hold;
    for n =1:sectionNumber
        subProfile = heartBoundary(n,:);
        L = sectionLength*Fs;
        Y = fft(subProfile);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        peaks = findpeaks(P1);
        beatFreqency(n) = f(find(P1==max(peaks)));
        normP1 = P1/max(peaks);
        
        %                  plot(f,normP1)
    end
    
%     mean(beatFreqency)
%     meanBeatFrequency(nn) = ans;
    % figure;
    % plot(f,P1)
    % title('Single-Sided Amplitude Spectrum of X(t)')
    % xlabel('f (Hz)')
    % ylabel('|P1(f)|')
    %      end
end

plot(beatFreqency);
clear dataCube;
save(strcat(folderName,'\process result'));