%determineShearAndScaling
%Use four lines on the edge of the projection to determine the Shear and
%Scaling for compensation.


% % Check Airy tails with:
% fileName='G:\2015-05-01 for de-con\cropped\16-07-15 beads 50ms stepSize162.5nm.avi';
% [pathstr,name,ext]=fileparts(fileName);
% if ~exist(strcat(pathstr,'\',name),'dir')
%     mkdir(fileName(1:end-4));
% end
% recStack=readDataCubeFromAviFile(fileName);
% imgSize=size(recStack);
% % functionName=mfilename();
% % configPath=mfilename('fullpath');
% % configPath=configPath(1:end-length(functionName));
% % defaultConfigFileName=fullfile(configPath,'OpenSPIM.json');
% % setupConfig=loadjson(defaultConfigFileName);
% dx=162.5e-9;
% xRange=dx*([1:imgSize(1)]-floor(imgSize(1)/2)-1);
% yRange=dx*([1:imgSize(2)]-floor(imgSize(2)/2)-1);
% zRange=dx*([1:imgSize(3)]-floor(imgSize(3)/2)-1);
% h=showImage(squeeze(max(recStack,[],1)).',-1,yRange*1e6,zRange*1e6);
% xlabel('horizontal [um]');ylabel('scan [um]')
% title('top projection');
% saveas(h,strcat(pathstr,'\',name,'\',name,' top projection.png'));
% figure
% showImage(squeeze(max(recStack,[],2)).',-1,xRange*1e6,zRange*1e6);
% h=xlabel('vertical [um]');ylabel('scan [um]')
% title('left projection');
% saveas(h,strcat(pathstr,'\',name,'\',name,' left projection.png'));


% 
function [scanShear scaling]=determineShearAndScaling(lines)

% store the start and end points (vertical, horizontal, scan) in two rows of a matrix
if nargin<1 || isempty(lines)
    lines(:,:,1)=[NaN -80.92 12.68; NaN -80.76 -3.087]; % Using NaN to indicate the unknown projection dimension coordinate
    lines(:,:,2)=[NaN 76.7 38.35; NaN 80.11 18.85]; % make it a 3D matrix
    lines(:,:,3)=[-74.42 NaN -18.69; -75.73 NaN -37.05];
    lines(:,:,4)=[75.56 NaN 32.5; 76.54 NaN 17.06];
    %         scanShear=[-0.026489 0.034449];
    %         scaling=[659.739403 718.135433];
    %         lines(:,:,1)=[NaN -54.39 -20.54; NaN -54.39 -28.49]; % Using NaN to indicate the unknown projection dimension coordinate
    %         lines(:,:,2)=[NaN 49.58 -0.74; NaN 50.5  -11.29]; % make it a 3D matrix
    %         lines(:,:,3)=[-40.33 NaN 24.97; -40.88 NaN 11.84];
    %         lines(:,:,4)=[38.85 NaN 18.13; 39.04 NaN 5.55];
    %         scanShear=[-0.012798 0.045418];
    %         scaling=[716.431715 835.045403];
    %         lines(:,:,1)=[NaN -58.46 -13.69; NaN -58.65 -27.20]; % Using NaN to indicate the unknown projection dimension coordinate
    %         lines(:,:,2)=[NaN 56.8 27.19; NaN 57.91 14.8]; % make it a 3D matrix
    %         lines(:,:,3)=[-41.99 NaN 5.18; -42.36 NaN -8.88];
    %         lines(:,:,4)=[43.29 NaN 1.85; 43.48 NaN -7.4];
    %         scanShear=[-0.003219 0.038299];
    %         scaling=[547.642941 894.245831];
    
    lines=lines*1e-6; % all in units of micron
end

isValid=~isnan(lines);
dataSize=size(lines);

positionsXY=mean(lines(:,1:2,:));
vectors=lines(2,:,:)-lines(1,:,:);
gradient=vectors(1,1:2,:)./repmat(vectors(1,3,:),[1 2 1]);

gradientNoNaN=gradient; gradientNoNaN(isnan(gradientNoNaN))=0;
positionsXYNoNaN=positionsXY; positionsXYNoNaN(isnan(positionsXYNoNaN))=0;
meanGradient=sum(gradientNoNaN,3)./sum(all(isValid(:,1:2,:)),3);
meanPositionsXY=sum(positionsXYNoNaN,3)./sum(all(isValid(:,1:2,:)),3);
relativePositionXY=(positionsXYNoNaN-repmat(meanPositionsXY,[1 1 dataSize(3)])).*all(isValid(:,1:2,:));

scaling=-sum((gradientNoNaN-repmat(meanGradient,[1 1 dataSize(3)])).*relativePositionXY,3)./sum(abs(relativePositionXY).^2,3);

scanShear=-meanGradient-scaling.*meanPositionsXY;

if nargout<1,
    logMessage('\nscanShear=[%0.6f %0.6f];\nscaling=[%0.6f %0.6f];',[scanShear, scaling]);
    clear scanShear;
end


end
