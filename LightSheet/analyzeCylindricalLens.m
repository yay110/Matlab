%
% Open Zemax before running this program
%
function analyzeCylindricalLens(lensFile)
    close all;
    if nargin<1
        lensFile='Z:/Tom/zemax/single lens airy light sheet - tom.zmx';
    end
    angles=[0:.5:60];
    testWavelengths=[405 488 532 561 633]*1e-9;
        
    zDDEInit();
    try
        zLoadFile(lensFile);
        zPushLens(1000,true);
        zSetWaveMatrix(cat(2,testWavelengths(:)*1e6,ones(size(testWavelengths(:)))));
        
        apertureDiameter=zGetSystemAper()*[0; 0; 1e-3];
        
        wavelengths=zGetWaveMatrix(); wavelengths=wavelengths(:,1)*1e-6;

        tilts=zeros(numel(angles),numel(wavelengths));
        defocuses=tilts;
        alphas=tilts;
        rmsResiduals=tilts;
        backFocalLengths=tilts;
        axisHeights=tilts;
        axisDir=tilts;
        for angleIdx=1:numel(angles),
            % Rotate cylindrical lens
            zSetSurfaceParameter(2,3,angles(angleIdx));
            for wavelengthIdx=numel(wavelengths):-1:1,
                % Switch wavelength
                zSetWaveMatrix([wavelengths(wavelengthIdx)*1e6,1]);
                % Refocus
                merit=zOptimize(100);
                zPushLens(1000,true); % Update the GUI so we can follow
                backFocalLengths(angleIdx,wavelengthIdx)=zGetSurfaceData(7,3)*1e-3; % record the back focal length
                % Find image position
                rayTraceData = zGetTrace(1,0,-1,0.0,1.0,0.0,0.0);
                axisHeights(angleIdx,wavelengthIdx)=rayTraceData(4)*1e-3;
                axisDir(angleIdx,wavelengthIdx)=rayTraceData(10);
                % Analyze wavefront
                [tangentialPupilPos, tangentialOPD, sagittalPupilPos, sagittalOPD]=getOPDData();
                [tilts(angleIdx,wavelengthIdx) defocuses(angleIdx,wavelengthIdx) alphas(angleIdx,wavelengthIdx) rmsResiduals(angleIdx,wavelengthIdx)]=analyzeWavefront(tangentialPupilPos,tangentialOPD);
            end
        end
        zSetWaveMatrix(cat(2,wavelengths(:)*1e6,ones(size(wavelengths(:)))));
                
        fig=figure('Position',[100 100 800 900]);
        ax(1,1)=subplot(3,2,1);
        lns=plot(angles,backFocalLengths*1000); setColors(lns,wavelengths);
        lab(1,1)=xlabel('lens angle (degree)'); lab(5,2)=ylabel('back focal length (mm)');
        ax(2,1)=subplot(3,2,2);
        lns=plot(angles,-axisHeights*1000); setColors(lns,wavelengths);
        lab(2,1)=xlabel('lens angle (degree)'); lab(6,2)=ylabel('position (mm)');
        ax(3)=subplot(3,2,3);
        lns=plot(angles,-alphas); setColors(lns,wavelengths);
        lab(3,1)=xlabel('lens angle (degree)'); lab(3,2)=ylabel('\alpha (\lambda)');
        ax(4)=subplot(3,2,4);
        lns=plot(angles,rmsResiduals); setColors(lns,wavelengths);
        lab(4,1)=xlabel('lens angle (degree)'); lab(4,2)=ylabel('residual (\lambda)');
        
        ax(5)=subplot(3,2,5);
        lns=plot(angles,tilts); setColors(lns,wavelengths);
        lab(5,1)=xlabel('lens angle (degree)'); lab(1,2)=ylabel('tilt (\lambda)'); 
        ax(6,1)=subplot(3,2,6);
        lns=plot(angles,defocuses); setColors(lns,wavelengths);
        lab(6,1)=xlabel('lens angle (degree)'); lab(2,2)=ylabel('defocus (\lambda)');
        
        set([ax; lab(:)],'FontSize',16,'FontWeight','bold');
        set(ax(1),'YLim',[0 50]); set(ax(2),'YLim',[0 5]);
        set(ax(3),'YLim',[-1 25]); set(ax(4),'YLim',[-1 25]);
        set(ax(5),'YLim',[-15 15]); set(ax(6),'YLim',[-5 5]);
        set(ax(:),'XLim',[min(angles) max(angles)],'XTick',[0:10:100]);
        set(ax,'LineWidth',2,'TickDir','out','Box','off');
        
        leg=legend(ax(3,1),cellfun(@(w) sprintf('%0.0f nm',w*1e9),num2cell(testWavelengths),'UniformOutput',false),'FontSize',12,'Location','NorthWest')
        
        drawnow();
        
        outputFileName=lensFile(1:end-4);
        saveas(fig,strcat(outputFileName,'.fig'));
        saveas(fig,strcat(outputFileName,'.png'));
        print(fig,strcat(outputFileName,'.eps'),'-depsc2');
        plot2svg(strcat(outputFileName,'.svg'),fig);
        
        save(strcat(outputFileName,'.mat'),'wavelengths','angles','axisHeights','backFocalLengths','tilts','defocuses','alphas','rmsResiduals');
    catch Exc
        % Make sure to always close the connection again
        zDDEClose();
        rethrow(Exc);
    end
    zDDEClose();
end

function setColors(lns,wavelengths)
    for lineIdx=1:numel(lns)
        set(lns(end-(lineIdx-1)),'Color',spectrumToRGB(wavelengths(lineIdx),.5,true));
    end
end

function [tilt defocus alpha rmsResidual]=analyzeWavefront(tangentialPupilPos,tangentialOPD)
    for wavelengthIdx=1:size(tangentialOPD,2),
        [fitobject,gof,output]=fit(tangentialPupilPos.',tangentialOPD(:,wavelengthIdx),'poly3');
        tilt(wavelengthIdx)=fitobject.p3;
        defocus(wavelengthIdx)=fitobject.p2;
        alpha(wavelengthIdx)=fitobject.p1;
        
        err=tangentialOPD(:,1)-fitobject(tangentialPupilPos.');
        rmsResidual(wavelengthIdx)=sqrt(mean(abs(err).^2));
    end
    
end

function [tangentialPupilPos, tangentialOPD, sagittalPupilPos, sagittalOPD, wavelengths]=getOPDData()
    folderName=fileparts(mfilename('fullpath'));
    tempFile=fullfile(folderName,'temp.txt');
    
    % Calc OPD
    zGetTextFile(tempFile,'Opd',fullfile(folderName,'OPD.CFG'), 1);
    fid=fopen(tempFile);
    % Find start of tangential data
    line=fgetl(fid); line=line(2:2:end);
    while ischar(line) && isempty(regexpi(line,'Pupil'))
        line = fgetl(fid); line=line(2:2:end);
    end
    wavelengths=regexp(line,'(\d+\.?\d+)','tokens');
    wavelengths=cellfun(@(e) str2double(e{1}),wavelengths)*1e-6;
    % Copy tangential data
    tangentialPupilPos=0;
    tangentialOPD=zeros(1,numel(wavelengths));
    line = fgetl(fid); line=line(2:2:end);
    rayIdx=1;
    while ischar(line) && isempty(regexpi(line,'field'))
        values=regexp(line,'(-?\d+\.?\d+)','tokens');
        if ~isempty(values),
            values=cellfun(@(e) str2double(e{1}),values);
            tangentialPupilPos(rayIdx)=values(1);
            tangentialOPD(rayIdx,:)=values(2:end);
            rayIdx=rayIdx+1;
        end
        line = fgetl(fid); line=line(2:2:end);
    end
    while ischar(line) && isempty(regexpi(line,'Pupil'))
        line = fgetl(fid); line=line(2:2:end);
    end
    % Read sagittal OPDs
    sagittalPupilPos=0;
    sagittalOPD=zeros(1,numel(wavelengths));
    line = fgetl(fid); line=line(2:2:end);
    rayIdx=1;
    while ischar(line)
        values=regexpi(line,'(-?\d+\.?\d+)','tokens');
        if ~isempty(values),
            values=cellfun(@(e) str2double(e{1}),values);
            sagittalPupilPos(rayIdx)=values(1);
            sagittalOPD(rayIdx,:)=values(2:end);
            rayIdx=rayIdx+1;
        end
        line = fgetl(fid); line=line(2:2:end);
    end        
    fclose(fid);
end