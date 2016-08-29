function dataCube = readImagesFromFolder(folderName,format)

% This function will read the images format defined in the folder and return a matrix cube
% example: dataCube = readImagesFromFolder(folderName,'.tiff')

if nargin < 2
    format = '.tiff';
end

if nargin < 1
    error('Folder that contains images should be defined!');
end

fileNames = dir(folderName);

cubeIdx = 1;

for n = 1:length(fileNames)
    if (~fileName(n).isdir) && strfind(fileName(n).name,format)
        dataCube(:,:,cubeIdx) = imread(strcat(folderName,'\',fileNames(n).name));
        cubeIdx = cubeIdx +1;
    end
end

end