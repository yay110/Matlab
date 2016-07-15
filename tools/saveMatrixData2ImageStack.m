function saveMatrixData2ImageStack(dataCube,outputFolder)
%save matrix in current workspace as .png files to be view
%in, e.g., Image J.  Input folders contain one or more
%deconvolved data cubes (the .mat file output of the processWaterImmersion
%algorithm) or subfolders containing them.

%Select data sets to convert.  Options are 'lightSheetPsf', 'recorded', and
%'deconvolved'
outputFolder = strcat(outputFolder,'_deconvolved');

logMessage('Saving deconvolved image stack to .png files....');

if ~exist(outputFolder,'dir');
    mkdir(outputFolder);
end


% 
% fileDir=dir(folderName);
% if length(fileDir)>2 %=2 means empty folder
%     for fileNo=3:length(fileDir)
%         fileName=fileDir(fileNo).name;
%         fullFileName=strcat(folderName,'\',fileName);
%         
%         %Regression to search nested folders for .mat files
%         if fileDir(fileNo).isdir==1 && strcmp(fileName,'ResultswithXoffset')==0
%             convertLSMatfileData2ImageStack(fullFileName,dataTypes);
%             
%             %Check current folder for .mat files
%         else
%             if strcmp(fileName(end-3:end),'.mat')
%                 
%                 %check matfile contains correct variables
%                 Flag_wrongMatfile=0;
%                 try
%                     load(fullFileName,'xRange');
%                     clear xRange
%                 catch err
%                     Flag_wrongMatfile=1;
%                 end
%                 
%                 %if correct matfile variables, load datacubes
%                 %and define folder(s) to save image stacks
%                 if ~Flag_wrongMatfile
%                     disp(strcat('Now converting file: ',fullFileName))
%                     
%                     %Handle data types one at a time
%                     for dataNo=1:length(dataTypes)
%                         if strcmp(dataTypes{dataNo},'lightSheetPsf')
%                             disp('Converting lightSheetPsf to image')
%                             load(fullFileName,'lightSheetPsf')
%                             dataCube=squeeze(lightSheetPsf)';
%                             clear lightSheetPsf
%                             outputFolder=strcat(fullFileName(1:end-4),'\lightSheetPsf');
%                         end
%                         if strcmp(dataTypes{dataNo},'deconvolved')
%                             disp('Converting deconvolved data to image stack')
%                             load(fullFileName,'restoredDataCube')
%                             dataCube=restoredDataCube;
%                             clear restoredDataCube
%                             outputFolder=strcat(fullFileName(1:end-4),'\deconvolvedImages');
%                         end
%                         if strcmp(dataTypes{dataNo},'recorded')
%                             disp('Converting recorded data to image stack')
%                             load(fullFileName,'recordedImageStack')
%                             dataCube=recordedImageStack;
%                             clear recordedImageStack
%                             outputFolder=strcat(fullFileName(1:end-4),'\recordedImages');
%                         end
%                         dataMax=max(dataCube(:));
%                         mkdir(outputFolder);


                         dataMax=max(dataCube(:));
                        %Write images of each designated
                        %data cube's XY slices to
                        %appropriate folder
                        for imageNo=1:size(dataCube,3)
                            outputFullFileName=strcat(outputFolder,'\','Image_',num2str(10000+imageNo),'.png');
                            outputImage=double(squeeze(dataCube(:,:,imageNo)));
%                             outputImage=outputImage.*uint16(outputImage>0);
                            outputImage=outputImage./double(dataMax);
                            imwrite(outputImage,outputFullFileName,'png','bitdepth',8);
                        end
%                         clear dataCube
                        logMessage('Deconvolved image stack written to file');
%                     end
%                 end
%             end
%         end
%     end
% end


end