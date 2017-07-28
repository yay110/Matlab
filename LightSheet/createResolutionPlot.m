function createResolutionPlot(inputFolder)
    close all;
    
    if nargin<1 || isempty(inputFolder)
        beamTypes={'Gaussian','Bessel10','Bessel5','Airy'};
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\forwardscans2';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\forwardscans';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\backwardscans';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\morePowerBackward';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\morePowerForward';

        inputFolder='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\resolutionVsFOV\UsedData';
    end
    
    requiredAspectRatioOfPlot=3.1678;
    
    outputFilePath='C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\resolutionVsFOV\tmp\';
    
    for beamTypeIdx=1:length(beamTypes),
        samples{beamTypeIdx}=load(fullfile(inputFolder,['images',beamTypes{beamTypeIdx},'.mat']));
        samples{beamTypeIdx}=samples{beamTypeIdx}.samples;
    end
    
    %Keep only those points that fall within the plot
    positions=cell2mat({samples{1}.centroid}.');
    inPlot=positions(:,1)>=-20e-6;
    for beamTypeIdx=1:length(beamTypes),
        samples{beamTypeIdx}=samples{beamTypeIdx}(inPlot);
    end
    positions=positions(inPlot,:);
    
    positions=cell2mat({samples{1}.centroid}.');
    for beamTypeIdx=1:length(beamTypes),
        fwhm{beamTypeIdx}=cell2mat({samples{beamTypeIdx}.fwhm}.');
    end
%     distanceFromWaist=-positions(:,2);
    distanceFromWaist=abs(positions(:,2));
    
    fig=figure('Position',[50 50 [800 300]*1.25],'NumberTitle','off','Name',outputFilePath);
    ax=axes();
    scatterMarkers{1}=scatterV6(distanceFromWaist*1e6,fwhm{1}(:,3)*1e6,40,[1 1 1],'Marker','+','MarkerEdgeColor',[0 0 0],'MarkerFaceColor','none','LineWidth',2);
    hold on;
    scatterMarkers{2}=scatterV6(distanceFromWaist*1e6,fwhm{2}(:,3)*1e6,50,[1 1 1],'Marker','x','MarkerEdgeColor',[0 0 0.75],'MarkerFaceColor','none','LineWidth',2);
    hold on;
    scatterMarkers{3}=scatterV6(distanceFromWaist*1e6,fwhm{3}(:,3)*1e6,40,[1 1 1],'Marker','^','MarkerEdgeColor',[1 0 0],'MarkerFaceColor','none','LineWidth',2);
    hold on;
    scatterMarkers{4}=scatterV6(distanceFromWaist*1e6,fwhm{4}(:,3)*1e6,40,[0 0.5 0],'filled','MarkerEdgeColor','none'); %,'filled'
    hold off;
    xlim([0 150]); ylim([0 25]);
    labs(1)=xlabel('x [\mum]'); labs(2)=ylabel('FWHM [\mum]');
    set(ax,'XTick',[0:25:1000]); set(ax,'YTick',[0:5:100]);
    set(ax,'LineWidth',2);
    set([ax labs],'FontSize',22,'FontWeight','bold','FontName','Arial');
    set(labs,'FontSize',24);
    set(ax,'Box','on');
    set(ax,'Units','pixel');
    axPos=get(ax,'Position');
    set(ax,'Position',[axPos(1:2) axPos(3)*[1 1/requiredAspectRatioOfPlot]]);
    drawnow();
    
    %attach callback
    allMarkers=cell2mat(scatterMarkers);
    set(allMarkers,'HitTest','on');
    for sampleIdx=1:size(allMarkers,1)
        for beamIdx=1:size(allMarkers,2);
            set(allMarkers(sampleIdx,beamIdx),'UserData',sampleIdx);
            sampleForEachBeam=cellfun(@(s)s(sampleIdx),samples);
            set(allMarkers(sampleIdx,beamIdx),'ButtonDownFcn',@(obj,evt) inspectSample(allMarkers(sampleIdx,:),beamTypes{beamIdx},sampleIdx,sampleForEachBeam));
        end
    end
    
    userData=struct();
    userData.outputFilePath=outputFilePath;
    userData.allMarkers=allMarkers;
    userData.samples=samples;
    userData.beamTypes=beamTypes;
    uicontrol(fig,'Style','pushbutton','String','Save','Callback',@saveFigure,'UserData',userData);
    
    set(fig,'ToolBar','figure');
    zoom(fig);
    
end

function markers=scatterV6(X,Y,S,varargin)
    nbMarkers=length(X);
    maxMarkers=50;
    markers=zeros(size(X));
    for markerIdx=1:maxMarkers:nbMarkers
        opts=varargin;
        selI=markerIdx:min(nbMarkers,markerIdx+maxMarkers-1);
        if all(ischar(S)) || isempty(opts) || all(isscalar(opts{1}))
            scatterGroup=scatter(X(selI),Y(selI),S);
        else
            if (length(opts)>1 && strcmpi(opts{2},'filled'))
                scatterGroup=scatter(X(selI),Y(selI),S,opts{1:2});
                opts={opts{3:end}};
            else
                scatterGroup=scatter(X(selI),Y(selI),S,opts{1});
                opts={opts{2:end}};
            end
        end
        for (argIdx=1:2:length(opts)-1)
            set(scatterGroup,opts{argIdx},opts{argIdx+1});
        end
        hold on;
        newMarkers=get(scatterGroup,'Children');
        markers(markerIdx-1+[1:length(newMarkers)])=newMarkers(end:-1:1); % Assume that the order is swapped!
    end
    hold off;
end

function saveFigure(obj,evt)
    fig=get(obj,'Parent');
    userData=get(obj,'UserData');
    
    % Save the data
    for markerIdx=1:size(userData.allMarkers,1)
        marker=userData.allMarkers(markerIdx,1);
        removed=strcmp(get(marker,'Visible'),'off');
        retained(get(marker,'UserData'))=~removed;
%         if removed
%             delete(userData.allMarkers(markerIdx,:));
%         end
    end
    for beamIdx=1:numel(userData.beamTypes)
        samples=userData.samples{beamIdx}(retained);
        save(fullfile(userData.outputFilePath,['images',userData.beamTypes{beamIdx},'.mat']),'samples');
    end
    
    % Save the figure
    description='';
    for beamTypeIdx=1:length(userData.beamTypes)
        description=[description,'_',userData.beamTypes{beamTypeIdx}];
    end
    fullFileNameFigure=fullfile(userData.outputFilePath,[description(2:end),'.svg']);
    
    set(obj,'Visible','off');
    logMessage('Writing figure to %s.',fullFileNameFigure);
    plot2svg(fullFileNameFigure,fig);
    set(obj,'Visible','on');
end

function inspectSample(sampleMarkers,sampleSet,sampleIdx,sampleForEachBeam)
    description=sprintf([sampleSet ' sample %d, centroid=(%0.1f,%0.1f,%0.1f)um, FWHM=%0.0f, %0.0f, %0.0f and %0.0f nm'],[sampleIdx,sampleForEachBeam(end).centroid*1e6 [sampleForEachBeam(1).fwhm(3) sampleForEachBeam(2).fwhm(3) sampleForEachBeam(3).fwhm(3) sampleForEachBeam(end).fwhm(3)]*1e9]);
    logMessage(description);
    
    colorMap=hot(1024);
    
    inspectFig=figure('Position',[200 200 [640 480]*1.5],'NumberTitle','off','Name','Check if this is a single bead.');
    for beamIdx=1:length(sampleForEachBeam)
        sample=sampleForEachBeam(beamIdx);
        subplot(length(sampleForEachBeam),3,1+3*(beamIdx-1));showImage(mapColor(sample.test.proj3./max(sample.test.proj3(:)),colorMap),[],sample.yRange*1e6,sample.xRange*1e6);axis equal;
        xlabel('x (propagation) [um]'); ylabel('y (swipe) [um]');
        subplot(length(sampleForEachBeam),3,2+3*(beamIdx-1));showImage(mapColor(sample.test.proj2./max(sample.test.proj2(:)),colorMap),[],sample.yRange*1e6,sample.zRange*1e6); axis equal;
        xlabel('x (propagation) [um]'); ylabel('z (scan) [um]');
        subplot(length(sampleForEachBeam),3,3+3*(beamIdx-1));showImage(mapColor(sample.test.proj1./max(sample.test.proj1(:)),colorMap),[],sample.xRange*1e6,sample.zRange*1e6); axis equal;
        xlabel('y (swipe) [um]'); ylabel('z (scan) [um]');
    end
    
    set(inspectFig,'NumberTitle','off','Name',description);
    set(inspectFig,'ToolBar','figure');
    zoom(inspectFig);
    
    drawnow();

    keepAction=@(obj,evt) close(inspectFig);
    
    uicontrol(inspectFig,'Style','pushbutton','String','Keep','Position',[10 10 50 20],'Callback',keepAction);
    uicontrol(inspectFig,'Style','pushbutton','String','Remove','Position',[80 10 50 20],'Callback',@(obj,evt) removeMarkers(sampleMarkers,inspectFig));
    
end
function removeMarkers(sampleMarkers,inspectFig)
    set(sampleMarkers,'Visible','off');
    close(inspectFig);
end
    