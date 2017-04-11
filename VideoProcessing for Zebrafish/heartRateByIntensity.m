%% =================process Images =================
clear;

dataFolder =["E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\1. before drug";...
    ''E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\2.adding and after drug'';...
    ''E:\2017-03-29 Zebrafish heartbeat\6th sample 2Vpp 1.457Mhz\3.washing away drug'';...
    ''E:\2017-03-29 Zebrafish heartbeat\7th sample\1. before drug'';...
    ''E:\2017-03-29 Zebrafish heartbeat\7th sample\2. adding to 320mgL'';...
    ''E:\2017-03-29 Zebrafish heartbeat\7th sample\3. adding to 500ml'';...
    ''E:\2017-03-29 Zebrafish heartbeat\7th sample\4. washing away'';...
    ''E:\2017-03-29 Zebrafish heartbeat\8th sample\1.before drug'';...
    ''E:\2017-03-29 Zebrafish heartbeat\8th sample\2. adding to 320mgL'';...
    ''E:\2017-03-29 Zebrafish heartbeat\8th sample\3. adding to 500mgL'';...
    ''E:\2017-03-29 Zebrafish heartbeat\8th sample\4. washing away''];

for folder = 1:%length(dataFolder)
    
    folderName = dataFolder(10)%folder);
    
    folderNames = dir(string(folderName));
    
    if folderNames(3).isdir
        for n = 1:length(folderNames)-2
            subFolderNames(n).name = strcat(folderName,'\',folderNames(n+2).name);
        end
        sectionLength = 9;
    else
        subFolderNames.name = folderName;
        sectionLength = 9;
    end
    
%     tic
    for nn = 1:length(subFolderNames)
        cd(subFolderNames(1).name)
        fileName = dir('**.tif*');
        
        Fs = 20;            %sampling frequency
        frameNumber = Fs*sectionLength;
        %read the images
        
        sectionNumber = floor(length(fileName)/(frameNumber));
        
        % read the first frame
        firstFrame = imread(strcat(subFolderNames(nn).name,'\',fileName(2).name));
        % find maxima on the frame
        [x y] = find(firstFrame ==max(firstFrame(:)));
        
        
        for n=1:sectionNumber
            for m = 1:frameNumber
                currentFrame = imread(strcat(subFolderNames(nn).name,'\',fileName((n-1)*frameNumber+m).name));
                profile(nn,n,m) = sum(sum(currentFrame(x-50:x,y-50:y)));
            end
        end
        
        
%         t = toc;
%         remainingTime = t/nn * (length(subFolderNames)-nn)/60
    end
    
    profile = squeeze(profile);
    
    save(strcat(folderName,'\extractedProfileFromMaxIntensityRegion'));
    
end

