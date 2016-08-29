folderName = 'E:\2016-08-13 Ciona\Ciona patch 3\stacks\1 stage';

directoryList = dir(folderName);

nrStacks = 1;
for listIdx = 1:length(directoryList)
    %     if (directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.' && strfind(directoryList(listIdx).name,'_stack'))
    if (directoryList(listIdx).isdir && ~isempty(strfind(directoryList(listIdx).name,'stack')))
        stackList(nrStacks) = directoryList(listIdx);
        nrStacks = nrStacks+1;
    end
end


for listIdx=60:length(stackList),
    expandedFolderName=strcat(folderName,'\',stackList(listIdx).name);
    fileNames = dir(expandedFolderName);
    
    cubeIdx = 1;
    for n = 1:length(fileNames)
        if (~fileNames(n).isdir) && strfind(fileNames(n).name,'.tiff')
            dataCube(:,:,cubeIdx) = imread(strcat(expandedFolderName,'\',fileNames(n).name));
            cubeIdx = cubeIdx + 1;
        end
    end
    for n = 1:length(dataCube(1,1,:))
        %         Zprofile(n) = sum(sum(dataCube(:,:,n)));
        %         Zmax(n) = max(max(dataCube(:,:,n)));
        conditionIdx(n) = sum(sum(dataCube(:,:,n) > 200));
    end
    %         %Normalize the cube Z profile to single pixel;
    %
    %         frameSize = size(dataCube);
    %         Zprofile = Zprofile/(frameSize(1)*frameSize(2));
    
    for n = 2:length(conditionIdx)
        if conditionIdx(n-1)<200 && conditionIdx(n)>200
            stackStartIdx = n;
            break;
        end
    end
    for n = stackStartIdx:length(conditionIdx)
        if conditionIdx(n-1)>200 && conditionIdx(n)<200
            stackEndIdx = n;
            break
        end
        
    end
    newdataCube = dataCube(:,:,stackStartIdx:stackEndIdx);
    
    maxProjection = max(newdataCube,[],3);
    imwrite(maxProjection,strcat(expandedFolderName,'.tiff'));
%     Zprojection(:,:,listIdx) = maxProjection;
%     magesc(double(maxProjection)/max(double(maxProjection(:))));
end