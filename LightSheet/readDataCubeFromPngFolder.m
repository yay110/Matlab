% [dataCube maxValue]=readDataCubeFromPngFolder(folderName,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     fileName: string representing the folder containing a series of png files numbered starting from 0.png
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%
% Outputs:
%     dataCube: the 3D matrix of values
%     maxValue: the maximum value that could have been stored in this file.
%
% Returns a 3D matrix in single precision and normalized to 1;
%
function [dataCube maxValue]=readDataCubeFromPngFolder(folderName,projectionDimensionOrSubCube,frameIndexes,normalize)
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    if (nargin<3)
        frameIndexes=-1;
    end
    frameIndexes=sort(frameIndexes);
    if (nargin<4 || isempty(normalize))
        normalize=true;
    end
    
    maxValue=2^8-1;
    getFrameFileName=@(frameIdx) fullfile(folderName,sprintf('%d.png',frameIdx));
    getFrame=@(frameIdx) mean(single(imread(getFrameFileName(frameIdx))),3);
    
    % Find out how many frames we have
    frameIdx=1;
    while exist(fullfile(folderName,sprintf('%d.png',frameIdx)),'file'),
        frameIdx=frameIdx+1;
    end
    nbFrames=frameIdx-1;
    
    firstFrame=getFrame(0);
    imgSize=size(firstFrame);
    
    if (length(frameIndexes)==1 && frameIndexes(1)==-1)
        frameIndexes=[1:nbFrames];
    else
        frameIndexes=intersect([1:nbFrames],frameIndexes);
    end
    
    dataCube=zeros([imgSize length(frameIndexes)],'single');
    for frameIndexIdx = 1:length(frameIndexes)
        frameIndex=frameIndexes(frameIndexIdx);
        if frameIndexIdx>1
            img=getFrame(frameIndex);
        else
            img=firstFrame;
        end
        % (project and) store
        if (isempty(projectionDimensionOrSubCube))
            %Complete data cube
            dataCube(:,:,frameIndexIdx)=img;
        else
            %Project or crop the data cube
            if (max(size(projectionDimensionOrSubCube))==1)
                %Project the full cube along one dimension specified by projectionDimensionOrSubCube
                if (any(projectionDimensionOrSubCube==1))
                    img=max(img,[],1);
                end
                if (any(projectionDimensionOrSubCube==2))
                    img=max(img,[],2);
                end
                if (~any(projectionDimensionOrSubCube==3))
                    dataCube(:,:,frameIndexIdx)=img;
                else
                    dataCube(:,:,1)=dataCube(:,:,1)+img;
                end
            else
                %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
                if(frameIdx>=projectionDimensionOrSubCube(3,1) && frameIdx<=projectionDimensionOrSubCube(3,2))
                     img=img(projectionDimensionOrSubCube(1,1):projectionDimensionOrSubCube(1,2),projectionDimensionOrSubCube(2,1):projectionDimensionOrSubCube(2,2));
                     dataCube(:,:,frameIndexIdx-projectionDimensionOrSubCube(3,1)+1)=img;
                end
            end
        end
    end
    if (normalize)
        dataCube=dataCube./maxValue;
    end
end