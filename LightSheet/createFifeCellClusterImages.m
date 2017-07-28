function createFifeCellClusterImages()
    close all;
    
    pathToData='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-2B';
    [projRawGaussianY projRawGaussianZ projDecGaussianY projDecGaussianZ xRange yRange zRange]=selectSubSetSingle([pathToData '/recording_lambda488nm_alpha0_beta100.mat']);
    [projRawAiryY projRawAiryZ projDecAiryY projDecAiryZ]=selectSubSetSingle([pathToData '/recording_lambda488nm_alpha2_beta100.mat']);
    
    outputFolder='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\FiveCellCluster';
    
    normalization=max([projRawGaussianY(:); projRawGaussianZ(:); projDecAiryY(:); projDecAiryZ(:)]);
%     projRawGaussianY=projRawGaussianY./normalization;
%     projRawGaussianZ=projRawGaussianZ./normalization;
%     projDecAiryY=projDecAiryY./normalization;
%     projDecAiryZ=projDecAiryZ./normalization;
    
    colorMap=@(values) interpolatedColorMap(values,[0 0 0; 0 1 0; 1 1 1],[0.06 .5 1]);
    projRawGaussianYColor=mapColor(projRawGaussianY./normalization,colorMap);
    projRawGaussianZColor=mapColor(projRawGaussianZ./normalization,colorMap);
    projDecAiryYColor=mapColor(projDecAiryY./normalization,colorMap);
    projDecAiryZColor=mapColor(projDecAiryZ./normalization,colorMap);
    
    fig=figure('Position',[50 50 900 900]);
    
    axs(1)=subplot(2,2,1);
    showImage(projRawGaussianZColor,[],yRange*1e6,xRange*1e6,axs(1));
    formatAxes(axs(1),'X (propagation) [\mum]','Y (swipe) [\mum]');
    axs(2)=subplot(2,2,2);
    showImage(projDecAiryZColor,[],yRange*1e6,xRange*1e6,axs(2));
    formatAxes(axs(2),'X (propagation) [\mum]','Y (swipe) [\mum]');
    axs(3)=subplot(2,2,3);
    showImage(projRawGaussianYColor,[],yRange*1e6,zRange*1e6,axs(3));
    formatAxes(axs(3),'X (propagation) [\mum]','Z (axial) [\mum]');
    axs(4)=subplot(2,2,4);
    showImage(projDecAiryYColor,[],yRange*1e6,zRange*1e6,axs(4));
    formatAxes(axs(4),'X (propagation) [\mum]','Z (axial) [\mum]');
    linkaxes(axs(1:2));
    linkaxes(axs(3:4));
    
    drawnow();
    
    logMessage('Saving in %s...',outputFolder);
    print(fig,'-dpsc2',[outputFolder, '/fiveCellClusterRaw.eps']);
    imwrite(projRawGaussianYColor,[outputFolder, '/projRawGaussianY.png']);
    imwrite(projRawGaussianZColor,[outputFolder, '/projRawGaussianZ.png']);
    imwrite(projDecAiryYColor,[outputFolder, '/projDecAiryY.png']);
    imwrite(projDecAiryZColor,[outputFolder, '/projDecAiryZ.png']);
    
    imwrite(repmat(permute(flipud(colorMap(1024)),[1 3 2]),[1 64]),[outputFolder, '/colorMap.png']);
end

function formatAxes(ax,xLabel,yLabel)
    axesColor=[0 0 0];
    axis(ax,'equal');
    set(ax,'TickDir','out','TickLength',[0.01 0.025]*3,'Color',axesColor,'XColor',axesColor,'YColor',axesColor,'LineWidth',3,'FontSize',14,'FontWeight','bold');
    set(ax,'XTick',[-1000:10:1000],'YTick',[-1000:10:1000]);
    xlabel('Parent',ax,xLabel,'LineWidth',3,'FontSize',16,'FontWeight','bold');
    ylabel('Parent',ax,yLabel,'LineWidth',3,'FontSize',16,'FontWeight','bold');
        
end

function [projRawY projRawZ projDecY projDecZ xRange yRange tRange]=selectSubSetSingle(inputFileName)
    proj=@(data,dim) max(data,[],dim);

    logMessage('Loading %s...',inputFileName);
    load(inputFileName);
    
%     ranges=[-17.5 42.5; -45 15; -25 35]*1e-6; % square
    ranges=[-20 40; -45 15; -25 35]*1e-6; % square, after center correction
    yRangeSel=yRange>=-ranges(1,2) & yRange<-ranges(1,1); % x and y swapped between matlab and paper!
    xRangeSel=xRange>=ranges(2,1) & xRange<ranges(2,2);
    zRangeSel=zRange>=ranges(3,1) & zRange<ranges(3,2);
    recordedImageStack=recordedImageStack(xRangeSel,yRangeSel,zRangeSel);
    restoredDataCube=restoredDataCube(xRangeSel,yRangeSel,zRangeSel);
    xRange=xRange(xRangeSel);
    yRange=yRange(yRangeSel);
    zRange=zRange(zRangeSel);
    tRange=tRange(zRangeSel);
    
%     background=0.0075;
%     recordedImageStack=recordedImageStack-background;
%     restoredDataCube=restoredDataCube-background;
    
    projRawY=squeeze(proj(recordedImageStack,1)).';
    projRawZ=proj(recordedImageStack,3);
    clear recordedImageStack;
    projDecY=squeeze(proj(restoredDataCube,1)).';
    projDecZ=proj(restoredDataCube,3);
end