
clear
%% =================process Images =================
dataFolder(1).name = 'E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\1. before drug';
dataFolder(2).name = 'E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\2.adding and after drug';
dataFolder(3).name = 'E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\3.washing away drug';
dataFolder(4).name = 'E:\2017-03-29 Zebrafish heartbeat\7th sample\1. before drug';
dataFolder(5).name = 'E:\2017-03-29 Zebrafish heartbeat\7th sample\2. adding to 320mgL';
dataFolder(6).name = 'E:\2017-03-29 Zebrafish heartbeat\7th sample\3. adding to 500ml';
dataFolder(7).name = 'E:\2017-03-29 Zebrafish heartbeat\7th sample\4. washing away';
dataFolder(8).name = 'E:\2017-03-29 Zebrafish heartbeat\8th sample\1.before drug';
dataFolder(9).name = 'E:\2017-03-29 Zebrafish heartbeat\8th sample\2. adding to 320mgL';
dataFolder(10).name = 'E:\2017-03-29 Zebrafish heartbeat\8th sample\3. adding to 500mgL';
dataFolder(11).name = 'E:\2017-03-29 Zebrafish heartbeat\8th sample\4. washing away';

for folder = 5:7%:length(dataFolder)
    clear profile;
    folderName = dataFolder(folder).name
    
    folderNames = dir(folderName);
    
    if folderNames(3).isdir
        for n = 1:length(folderNames)-2
            subFolderNames(n).name = strcat(folderName,'\',folderNames(n+2).name);
        end
        sectionLength = 9;
    else
        subFolderNames.name = folderName;
        sectionLength = 60;
    end
    
    Fs = 20;            %sampling frequency
    frameNumber = Fs*sectionLength;
%     profile = zeros(length(subFolderNames),sectionNumber,frameNumber);
    
    for nn = 1:length(subFolderNames) 
        
        cd(subFolderNames(1).name)
        fileName = dir('**.tif*');        
        sectionNumber = ceil(length(fileName)/(frameNumber));

        %read the images          
        for n=1:sectionNumber
            n/sectionNumber
            % read the first frame
            firstFrame = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+1).name));
            % find maxima on the frame
            [x y] = find(firstFrame ==max(firstFrame(:)));
            
            for m = 1:frameNumber
                try
                currentFrame = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+m).name));
                catch
                    break
                end
                profile(nn,n,m) = sum(sum(currentFrame(x-50:x,y-50:y)));
                %                 profile1(nn,n,m) = sum(sum(currentFrame(x-50:x,y:y+50)));
                %                 profile2(nn,n,m) = sum(sum(currentFrame(x:x+50,y-50:y)));
                %                 profile3(nn,n,m) = sum(sum(currentFrame(x:x+50,y:y+50)));
                %                 [x y] = find(currentFrame ==max(currentFrame(:)));
                %                 profile4(nn,n,m) = x(1);
                %                 profile5(nn,n,m) = y(1);
            end            
        end        
    end
    
    profile = squeeze(profile);
    %     profile1 = squeeze(profile1);
    %     profile2 = squeeze(profile2);
    %     profile3 = squeeze(profile3);
    %     profile4 = squeeze(profile4);
    %     profile5 = squeeze(profile5);
    
    save(strcat(folderName,'\extractedProfileFromMaxIntensityRegion'));
    
end

