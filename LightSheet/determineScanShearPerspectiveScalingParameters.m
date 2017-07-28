function [dataStruct]=determineScanShearPerspectiveScalingParameters(recordedImageStack,xRange,yRange,zRange)
    % From a 3d datacube "recordedImageStack", allows determination of
    % scanShear and perspectiveScaling parameters by plotting the angle of
    % Airy PSF tails relative to the z-axis as a function of x and y.
    % Use on datasets of fluorescent beads only.
    
    %Manual file input in the case that the target .mat file is not in the
    %workspace.
    if nargin<1
        fileName='E:\LSM Data\2015-09-29_BeadTest\Aperture=0.6\2015-09-29 09_50_50.507\recording0_lambda532nm_alpha7_beta100.mat';
        load(fileName,'recordedImageStack','xRange','yRange','zRange');
    end
    
    %Unit conversion from µm to m
    xRange=xRange*1e6;
    yRange=yRange*1e6;
    zRange=zRange*1e6;
    
    xSpacing=xRange(2)-xRange(1);
    ySpacing=yRange(2)-yRange(1);
    zSpacing=zRange(2)-zRange(1);
    
    
    % Generate initial grid for localised searching for beads.
    % Let x and y refer to MATLAB coordinate system for clarity.
    noDataPoints=2000; %Number of random samples across image
    xStartingCoords=rand([1,noDataPoints],'single');
    yStartingCoords=rand([1,noDataPoints],'single');
    zStartingCoords=rand([1,noDataPoints],'single');
    xStartingCoords=round(1+(length(xRange)-1)*xStartingCoords);
    yStartingCoords=round(1+(length(yRange)-1)*yStartingCoords);
    zStartingCoords=round(1+(length(zRange)-1)*zStartingCoords);
    
    
    % Generate a data structure for storing data and populate some fields;
    dataStruct=[];
    for n=1:noDataPoints;
       dataStruct(n).beadNumber=n;
       dataStruct(n).startingXCoord=xStartingCoords(n);
       dataStruct(n).startingYCoord=yStartingCoords(n);
       dataStruct(n).startingZCoord=zStartingCoords(n);
       dataStruct(n).searchRangeX=[max(1,dataStruct(n).startingXCoord-50),min(length(xRange),dataStruct(n).startingXCoord+50)];
       dataStruct(n).searchRangeY=[max(1,dataStruct(n).startingYCoord-50),min(length(yRange),dataStruct(n).startingYCoord+50)];
       dataStruct(n).searchRangeZ=[max(1,dataStruct(n).startingZCoord-50),min(length(zRange),dataStruct(n).startingZCoord+50)];
       dataStruct(n).Flag_noBead=0;
       dataStruct(n).maxXCoord=[];
       dataStruct(n).maxYCoord=[];
       dataStruct(n).maxZCoord=[];
       dataStruct(n).analysisRangeX=[];
       dataStruct(n).analysisRangeY=[];
       dataStruct(n).analysisRangeZ=[];
       dataStruct(n).Flag_duplicatedBead=0;
       dataStruct(n).psfXZTilt=[];
       dataStruct(n).psfYZTilt=[];
    end
    clear xStartingCoords yStartingCoords zStartingCoords
    
%Randomly sample the image using the starting coordinates defined above
%(lines 25-27), searching for local maxima corresponding to the centers of
%bead PSFs.
    for n=1:noDataPoints;
        croppedDataCube=recordedImageStack(dataStruct(n).searchRangeX(1):dataStruct(n).searchRangeX(2),dataStruct(n).searchRangeY(1):dataStruct(n).searchRangeY(2),dataStruct(n).searchRangeZ(1):dataStruct(n).searchRangeZ(2));
        meanVal=mean(croppedDataCube(:));
        
        [maxVal,maxPos]=max(croppedDataCube(:));
        if maxVal==0;
            dataStruct(n).Flag_noBead=1;
        elseif maxVal<2*meanVal;
            dataStruct(n).Flag_noBead=1;
        else
            dataStruct(n).maxZCoord=floor(maxPos/size(croppedDataCube,2)/size(croppedDataCube,1))+dataStruct(n).searchRangeZ(1);
            if rem(maxPos,size(croppedDataCube,2)*size(croppedDataCube,1))==0;
                maxPos=size(croppedDataCube,2)*size(croppedDataCube,1);
                dataStruct(n).maxZCoord=dataStruct(n).maxZCoord-1;
            else
                maxPos=rem(maxPos,size(croppedDataCube,2)*size(croppedDataCube,1));
            end
            dataStruct(n).maxYCoord=floor(maxPos/size(croppedDataCube,1))+dataStruct(n).searchRangeY(1);
            if rem(maxPos,size(croppedDataCube,1))==0;
                maxPos=size(croppedDataCube,1);
                dataStruct(n).maxYCoord=dataStruct(n).maxYCoord-1;
            else
                maxPos=rem(maxPos,size(croppedDataCube,1));
            end
            dataStruct(n).maxXCoord=maxPos-1+dataStruct(n).searchRangeX(1);
            
%Define the analysis range, or volume in pixels of the region centered on
%the max coordinate that will be searched to define and fit the Airy tail.
%The analysis range can be changed to reflect the size of the bead PSFs at
%the magnification of the system.
            dataStruct(n).analysisRangeX=[max(1,dataStruct(n).maxXCoord-10),min(length(xRange),dataStruct(n).maxXCoord+10)];
            dataStruct(n).analysisRangeY=[max(1,dataStruct(n).maxYCoord-10),min(length(yRange),dataStruct(n).maxYCoord+10)];
            dataStruct(n).analysisRangeZ=[max(1,dataStruct(n).maxZCoord-20),min(length(zRange),dataStruct(n).maxZCoord+5)];
        end
    end
    
    for n=1:noDataPoints;
        if dataStruct(n).Flag_noBead~=1;
            % Plot CoM and Max value in each z-plane to determine PSF tilt     
            croppedDataCube=recordedImageStack(dataStruct(n).analysisRangeX(1):dataStruct(n).analysisRangeX(2),dataStruct(n).analysisRangeY(1):dataStruct(n).analysisRangeY(2),dataStruct(n).analysisRangeZ(1):dataStruct(n).analysisRangeZ(2));
            zVector=[1:size(croppedDataCube,3)];
            for IdZ=1:size(croppedDataCube,3)
                imagePlane=croppedDataCube(:,:,IdZ);
                meanVal=mean(imagePlane(:));
                [maxVal,maxPos]=max(imagePlane(:));
                if maxVal<2*meanVal
                    maxYPos(IdZ)=NaN;
                    maxXPos(IdZ)=NaN;
                    zVector(IdZ)=NaN;
                else
                    maxYPos(IdZ)=floor(maxPos/size(imagePlane,1))+1;
                    if rem(maxPos,size(imagePlane,1))==0
                        maxPos=size(croppedDataCube,1);
                        maxYPos(IdZ)=maxYPos(IdZ)-1;
                    else
                        maxPos=rem(maxPos,size(imagePlane,1));
                    end
                    maxXPos(IdZ)=maxPos;             
                end
                guessImage=zeros(size(imagePlane'));
                if isnan(maxXPos(IdZ))==0
                    guessImage(maxXPos(IdZ),maxYPos(IdZ))=1;
                end
            end
            
            for IdZ=1:size(croppedDataCube,3);
                if isnan(maxXPos(IdZ))
                    maxXPos(IdZ)=0;
                    maxYPos(IdZ)=0;
                    zVector(IdZ)=0;
                end
            end
            maxXPos=nonzeros(maxXPos);
            maxYPos=nonzeros(maxYPos);
            zVector=nonzeros(zVector);
            
            %Perform linear fitting to determine tilt
            pX=polyfit(zVector*zSpacing,maxXPos*xSpacing,1);
            pY=polyfit(zVector*zSpacing,maxYPos*ySpacing,1);
            fitX=polyval(pX,zVector*zSpacing);
            fitY=polyval(pY,zVector*zSpacing);
            dataStruct(n).psfXZTilt=pX(1);
            dataStruct(n).psfYZTilt=pY(1);
        
        elseif dataStruct(n).Flag_noBead==1;
            dataStruct(n).psfXZTilt=NaN;
            dataStruct(n).psfYZTilt=NaN;
        end
        clear zVector maxXPos maxYPos
    end

    %Store results
    for m=1:length(dataStruct);
        if dataStruct(m).Flag_noBead==0;
            xPosition(m)=dataStruct(m).maxXCoord;
            yPosition(m)=dataStruct(m).maxYCoord;
            xGradient(m)=dataStruct(m).psfXZTilt;
            yGradient(m)=dataStruct(m).psfYZTilt;
        else
            xPosition(m)=NaN;
            yPosition(m)=NaN;
            xGradient(m)=NaN;
            yGradient(m)=NaN;
        end
    end

    xPosition=xPosition(isnan(xPosition)~=1);
    yPosition=yPosition(isnan(yPosition)~=1);
    xGradient=xGradient(isnan(xGradient)~=1);
    yGradient=yGradient(isnan(yGradient)~=1);
    
    %Perform linear fit on distribution of tilt values to determine
    %correction parameters
    pXGrad=polyfit((xPosition-floor(length(xRange)/2))*xSpacing,xGradient,1);
    pYGrad=polyfit((yPosition-floor(length(yRange)/2))*ySpacing,yGradient,1);
    XGradfit=polyval(pXGrad,(xPosition-floor(length(xRange)/2))*xSpacing);
    YGradfit=polyval(pYGrad,(yPosition-floor(length(yRange)/2))*ySpacing);
    
    %Define output values
    scanShear=[pYGrad(2) pXGrad(2)]*-1
    scaling=[pYGrad(1) pXGrad(1)]*-1
    
    %Plot distributions of gradients (in XZ and YZ) and their linear fits
    figure(5);scatter((xPosition-floor(length(xRange)/2))*xSpacing,xGradient);
    xlabel('x [um]');ylabel('tilt angle [rad]');
    hold on;
    plot((xPosition-floor(length(xRange)/2))*xSpacing,XGradfit,'r');
    hold off;
    figure(6);scatter((yPosition-floor(length(yRange)/2))*ySpacing,yGradient);
    xlabel('y [um]');ylabel('tilt angle [rad]');
    hold on;
    plot((yPosition-floor(length(yRange)/2))*ySpacing,YGradfit,'r');
    hold off;
    
end