
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

for folder = 4:length(dataFolder)
    clear profile;
    folderName = dataFolder(folder).name%folder);
    
    folderNames = dir(folderName);
    
    if folderNames(3).isdir
        for n = 1:length(folderNames)-4
            subFolderNames(n).name = strcat(folderName,'\',folderNames(n+2).name);
        end
        sectionLength = 9;
    else
        subFolderNames.name = folderName;
        sectionLength = 9;
    end
    
    for nn = 1:length(subFolderNames)
        cd(subFolderNames(1).name)
        fileName = dir('**.tif*');
        
        Fs = 20;            %sampling frequency
        frameNumber = Fs*sectionLength;
        %read the images
        
        sectionNumber = floor(length(fileName)/(frameNumber));
        
        for n=1:sectionNumber
            
            % read the first frame
            firstFrame = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+1).name));
            % find maxima on the frame
            [x y] = find(firstFrame ==max(firstFrame(:)));
            
            for m = 1:frameNumber
                currentFrame = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+m).name));
                profile(nn,n,m) = sum(sum(currentFrame(x-50:x,y-50:y)));
            end
            
        end
        
    end
    
    profile = squeeze(profile);
    
    save(strcat(folderName,'\extractedProfileFromMaxIntensityRegion'));
    
end

