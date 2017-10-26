classdef Andor < handle
    % Andor class
    %
    % constructor Andor(CCDSettingTemp, AcquisitionMode, ExposureTime, ReadMode, TriggerMode, PreAmpGain)
    %               CCDSettingTemp, by default (-45 Celsius)
    %               AcquisitionMode, 1 - Single Scan, 2 - Accumulate, 3 - Kinetics, 4 - Fast Kinetics, 5 - Run till abort;
    %               ExposureTime, in seconds;
    %               ReadMode, 0 - Full Vertical Binning, 1 - Multi-Track, 2  - Random-Track, 3 - Single-Track, 4 - Image;
    %               TriggerMode,
    %         0 - internal, by defaulty
    %         1 - External
    %         6 - External Start
    %         7 - External Exposure (Bulb)
    %         9 - External FVB EM (only valid for EM Newton models in FVB mode)	10.	Software Trigger
    %         12 - External Charge Shifting
    %               PreAmpGain, 0 - 1x, 1 - 2x, 2 - 4x;
    %               NumberKinetics;             number of kinectcs; = 1 by default;
    %               NumberAccumulations;        Number of Accumulations = 1 by default;
    %               CCD_Pixels;                 %[CCD_XPixels, CCD_YPixels];   physical CCD pixels
    %               cropHeight;         CropMode, by default no crop, cropHeight = 0;
    %               CosmicRayFilterOn;      Cosmic Ray removal; 1 - on; 0 - off
    %
    %               CCDCurrentTemp;  %current CCD temperature;
    %               noGains;         %Available amplify gains
    %               CCDsize = []; % [CCD_XPixels CCD_YPixels];    %physical CCD pixels;
    %Shamrock
    %               ShamrockDev;     %Shamrock devices;
    %               Grating;         %grating used;
    %               SlitWidth;       %by default 150um;
    %               Wavelength;      %central wavelength;
    %               AxisWavelength   %Calibrated axis in nm;
    %
    %
    %
    % by default, andor();
    %
    % Example
    %       Andor(CCDSettingTemp, AcquisitionMode, ExposureTime, ReadMode, TriggerMode, PreAmpGain, NumberKinetics, NumberAccumulations,CropHeight,CosmicRayFilterOn)
    %
    %
    % Copyright @ Mingzhou Chen, @University of St. Andrews. ingzhou.chen@st-andrews.ac.uk, June 2016, Ver. 1.01
    properties
        CCDSettingTemp;
        AcquisitionMode;
        ExposureTime;
        ReadMode;
        TriggerMode;
        PreAmpGain;
        NumberAccumulations;
        NumberKinetics;
        CropHeight;   %CCD crop mode
        CosmicRayFilterOn; %Cosmic Ray filter on = 1 by default; 0 - off;
        
        SlitWidth;          %by default 150um;
        CurrentGrating;     %current grating in use;
        CentralWavelength;  %central wavelength;
        
        AxisWavelength   %Calibrated axis in nm;
        noGratings;         %Available gratings;
        
        CCDCurrentTemp;  %current CCD temperature;
        noGains;         %Available amplify gains
        
        CCD_Pixels = []; % [CCD_XPixels CCD_YPixels];    %CCD pixels in use;
    end
    properties(Dependent = true)
    end
    properties (SetAccess = protected)
    end
    properties (Access = protected)
        CCDsize = [];    %[xpixels ypixels]; %Physical CCD pixels size;
        ShamrockDev;     %Shamrock devices; only one shamrock dev can be controlled now
    end
    
    methods
        function Andor = Andor()
            fprintf('Initialising Andor camera... please wait!......');
            [ret]=AndorInitialize('');      %initialization camera
            CheckError(ret);
            if ret == atmcd.DRV_SUCCESS
                fprintf('done\nAndor camera has been initialized successfully!\n');
            else
                disp('Error occurs during Andor initialization! Please check...!')
            end
            
            [ret, iCoolerStatus] = IsCoolerOn();
            CheckWarning(ret);
            if ~iCoolerStatus
                [ret]=CoolerON();                             %   Turn on temperature cooler
                CheckWarning(ret);
                SetTemperature(-45);
                fprintf('Andor cooller is turned on and set to -45 degrees!\n');
            end
            
            %
            %             % Initialize parameters;
            %             [ret, mintemp, maxtemp] = GetTemperatureRange();
            %             CheckWarning(ret);
            %             if nargin < 1
            %                 CCDSettingTemp = -45;   %CCD temperature at -45 Celsius by default;
            %             elseif CCDSettingTemp<mintemp || CCDSettingTemp>maxtemp
            %                 CCDSettingTemp = -45;
            %             end
            %             Andor.CCDSettingTemp = CCDSettingTemp;
            %
            %get physical parameters of system;
            Andor.ReadMode = 4; %ReadMode, 0 - Full Vertical Binning, 1 - Multi-Track, 2  - Random-Track, 3 - Single-Track, 4 - Image;
%             [ret] = SetCropMode(0, 1, 0);
%             CheckWarning(ret);
            [ret,XPixels, YPixels]=GetDetector;
            CheckWarning(ret);
             Andor.CCDsize = [XPixels, YPixels];
            [ret]=SetImage(1, 1, 1, XPixels, 1, YPixels); %   Set the image size
            CheckWarning(ret);
            
            Andor.AcquisitionMode = 1;%AcquisitionMode,1 - Single Scan by default, 2 - Accumulate, 3 - Kinetics, 4 - Fast Kinetics, 5 - Run till abort;
            
            Andor.ExposureTime = 0.1;%   ExposureTime in second;
                        
            Andor.TriggerMode = 0;    %internal trigger mode
            
%             if ret == 20002
%                 if nargin < 6
%                     PreAmpGain = 0;                   % PreAmpGain, 0 - 1x by default, 1 - 2x, 2 - 4x;
%                 end
%                 Andor.PreAmpGain = PreAmpGain;
%             end
            [ret]=SetShutter(1, 1, 0, 0); % open shutter;
            CheckWarning(ret);
%             
%             if nargin < 7 %NumberKinetics
%                 NumberKinetics = 1; %by default
%             end
            Andor.NumberKinetics = 10;
%             
%             if nargin <8 %NumberAccumulations
%                 NumberAccumulations = 1; %by default
%             end
%             Andor.NumberAccumulations = NumberAccumulations;
                           
            fprintf('Andor is ready to use...\n');
        end
        
        function Andor = set.NumberAccumulations(Andor,newNumberAccumulations)
            if newNumberAccumulations>0
                [ret]=SetNumberAccumulations(newNumberAccumulations);
                CheckWarning(ret);
                Andor.NumberAccumulations = newNumberAccumulations;
            end
        end
        function Andor = set.NumberKinetics(Andor,newNumberKinetics)
            if newNumberKinetics>0
                [ret]=SetNumberKinetics(newNumberKinetics);
                CheckWarning(ret);
                Andor.NumberKinetics = newNumberKinetics;
            end
        end
        function Andor = set.CCDSettingTemp(Andor,newCCDSettingTemp)  %Set CCD temperature;
            [ret, CurrentTemp] = GetTemperature(); %get the CCD tempertature;
            if CurrentTemp == newCCDSettingTemp && ret == atmcd.DRV_TEMP_STABILIZED
                return;
            end
            [ret, iCoolerStatus] = IsCoolerOn();
            CheckWarning(ret);
            if ~iCoolerStatus
                [ret]=CoolerON();                             %   Turn on temperature cooler
                CheckWarning(ret);
            end
            
            [ret, mintemp, maxtemp] = GetTemperatureRange();
            CheckWarning(ret);
            if newCCDSettingTemp<mintemp || newCCDSettingTemp>maxtemp
                disp('Temperature is out of range! Set to -45 bt default!');
                newCCDSettingTemp = -45;
            end
            
            [ret, temp] = GetTemperature(); %get the CCD tempertature;
            if newCCDSettingTemp ~= temp || ret ~= atmcd.DRV_TEMPERATURE_STABILIZED
                ret = SetTemperature(newCCDSettingTemp);
                CheckWarning(ret);
                %waiting for cooling down and temp stablization;
                fprintf('Cooling Camera ---> ');
                msg = [num2str(temp) ' Celsius......'];
                mnum = length(msg);
                fprintf(msg);
                [ret, temp] = GetTemperature();
                while ret ~= atmcd.DRV_TEMPERATURE_STABILIZED % temperature is not stablized. Wait until set;
                    for mm = 1:mnum
                        fprintf('\b');
                    end
                    msg = [num2str(temp) ' Celsius......'];
                    mnum = length(msg);
                    fprintf(msg);
                    pause(1);
                    [ret, temp] = GetTemperature();
                    if ret == atmcd.DRV_NOT_INITIALIZED || ret == atmcd.DRV_ACQUIRING || ret == atmcd.DRV_ERROR_ACK || ret == atmcd.DRV_TEMPERATURE_OFF
                        %                         CheckWarning(ret);
                        warning('Temperature seeting error occurs...Please check!')
                        break;
                    end
                end
                Andor.CCDSettingTemp = newCCDSettingTemp;
                fprintf('...done.\n');
            end
        end
        function CCDCurrentTemp = get.CCDCurrentTemp(Andor)
            [ret, CCDCurrentTemp] = GetTemperature(); %get the CCD tempertature;
        end
        function Andor = set.AcquisitionMode(Andor,AcquisitionMode)
            AcquisitionMode = round(AcquisitionMode);
            if AcquisitionMode<1 || AcquisitionMode > 5 %AcquisitionMode,1 - Single Scan by default, 2 - Accumulate, 3 - Kinetics, 4 - Fast Kinetics, 5 - Run till abort;
                AcquisitionMode = 1;
            end
            [ret]=SetAcquisitionMode(AcquisitionMode);    %   Set acquisition mode; 1 for Single Scan
            CheckWarning(ret);
            Andor.AcquisitionMode = AcquisitionMode;
            switch AcquisitionMode
                case 1
                    disp('Using Single Scan mode now');
                case 2
                    disp('Using Accumulate mode now');
                case 3
                    disp('Using Kinetics mode now');
                case 4
                    disp('Using Fast Kinetics mode now');
                case 5
                    disp('Using Run till abort mode now');
                otherwise
            end
        end
        function Andor = set.ExposureTime(Andor,newExposureTime)
            [ret, MaxExp] = GetMaximumExposure();
            CheckWarning(ret);
            if newExposureTime>MaxExp
                newExposureTime = MaxExp;
            elseif newExposureTime<0.00090
                newExposureTime = 0.0009; %0.87mS in minimum exposure;
            end
            [ret]=SetExposureTime(newExposureTime); %   Set exposure time in seconds
            CheckWarning(ret);
            Andor.ExposureTime = newExposureTime;
        end
        function Andor = set.ReadMode(Andor,ReadMode)
            ReadMode = round(ReadMode);
            if ReadMode<0 || ReadMode > 4 %ReadMode, 0 - Full Vertical Binning, 1 - Multi-Track, 2  - Random-Track, 3 - Single-Track, 4 - Image;
                ReadMode = 0;
            end
            [ret]=SetReadMode(ReadMode);          %   Set read mode;
            CheckWarning(ret);
            Andor.ReadMode = ReadMode;
            switch ReadMode
                case 0
                    disp('Using FVB mode now');
                case 1
                    disp('Using Multi-Track mode now');
                case 2
                    disp('Using Random-Track mode now');
                case 3
                    disp('Using Single-Track mode now');
                case 4
                    disp('Using Image mode now');
                otherwise
            end
        end
        function Andor = set.TriggerMode(Andor,TriggerMode)
            if TriggerMode<0 || TriggerMode > 12
                %TriggerMode,
                %         0 - internal, by defaulty
                %         1 - External
                %         6 - External Start
                %         7 - External Exposure (Bulb)
                %         9 - External FVB EM (only valid for EM Newton models in FVB mode)	10.	Software Trigger
                %         12 - External Charge Shifting
                TriggerMode = 0;
            end
            [ret]=SetTriggerMode(TriggerMode);                   %   Set internal trigger mode by default
            if ret ~= 20002
                CheckWarning(ret);
                disp('Set to internal trigger by default');
                [ret]=SetTriggerMode(TriggerMode);
                CheckWarning(ret);
            end
            Andor.TriggerMode = TriggerMode;
        end
        function Andor = set.PreAmpGain(Andor,newPreAmpGain)
            if newPreAmpGain<0 || newPreAmpGain > 2 % PreAmpGain, 0 - 1x by default, 1 - 2x, 2 - 4x;
                newPreAmpGain = 0;
            end
            [ret]=SetPreAmpGain(newPreAmpGain);          %   Set read mode;
            CheckWarning(ret);
            Andor.PreAmpGain = newPreAmpGain;
        end
        function CCD_Pixels = get.CCD_Pixels(Andor)
            [ret,XPixels, YPixels]=GetDetector();           %   Get the CCD size
            CheckWarning(ret);
            CCD_Pixels = [XPixels YPixels];
        end
        function Andor = set.CropHeight(Andor,newCropHeight)
            Andor.CropHeight = setCropParam(Andor,newCropHeight);
        end
        function noGains = get.noGains(Andor)
            [ret, noGains] = GetNumberPreAmpGains();
            CheckWarning(ret);
        end
        function Andor = set.CosmicRayFilterOn(Andor,newCosmicRayFilterOn)
            if newCosmicRayFilterOn>0
                Andor.CosmicRayFilterOn=1;
                [ret] = SetFilterMode(2);
                CheckWarning(ret);
            else
                Andor.CosmicRayFilterOn=0;
                [ret] = SetFilterMode(0);
                CheckWarning(ret);
            end
        end
        function CosmicRayFilterOn = get.CosmicRayFilterOn(Andor)
            [ret,mode] = GetFilterMode();
            CheckWarning(ret);
            if mode>0
                CosmicRayFilterOn=1;
            else
                CosmicRayFilterOn=0;
            end
        end
    end
    
    methods
        function saveSIF(Andor,file,comment)
            [ret,gstatus]=AndorGetStatus();
            CheckWarning(ret);%
            if(gstatus ~= atmcd.DRV_IDLE)
                disp('System is busy, please try it again later.');
                return;
            end
            if nargin>2
                [ret] = SetSifComment(comment);
                CheckWarning(ret);
            end
            ret = SaveAsSif(file);
            CheckWarning(ret);
            
            %strange! SaveAsSif function because it save a file named C or D into the current folder; have to move it to the distination.
            movefile([pwd filesep file(1)],file); %alwasy save the filename as the first letter; need to move it to the distination;
            
            x = 1:Andor.CCD_Pixels(1);
            B = [ones(1,Andor.CCD_Pixels(1));x;x.^2;x.^3].';
            wavelengthInfo = pinv(B)*Andor.AxisWavelength;
            newLine = sprintf('%4.12f %3.12f %3.12e %3.12e',wavelengthInfo);
            replaceWavelengthInfoInSif(file,newLine); %add wavelength informatin into sif;
        end
        function acquireLive(Andor)   %
            fprintf('Starting live acquisition (Close the figure or press 0 to exit).......\n');
            
            YPixels = Andor.CCD_Pixels(2);
            XPixels = Andor.CCD_Pixels(1);
            if Andor.ReadMode == 4 %image mode
                I = zeros(YPixels,XPixels);
                fig = figure(1000);
                h = imagesc(I);axis image;
                colormap(gray);
            elseif Andor.ReadMode == 0 %FAB mode;
                I = zeros(XPixels,1);
                fig = figure(1000);
                h = plot(Andor.AxisWavelength, I);
            end
            [ret]=SetShutter(1, 1, 0, 0); %close shutter
            CheckWarning(ret);
            
            [ret] = PrepareAcquisition();
            dcm_obj = datacursormode(fig);
            set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
            set(dcm_obj,'UpdateFcn',@myupdatefcn);
            c_info = [];
            old_info = [];
            
            key = get(fig,'CurrentCharacter');
            set(fig,'Name','Real-time acquiring now......Close the figure or press 0 to exit');
            while not(strcmp(key,'0'))
                [ret] = StartAcquisition();
                CheckWarning(ret);
                [ret,gstatus]=AndorGetStatus;
                CheckWarning(ret);
                while(gstatus ~= atmcd.DRV_IDLE)
                    pause(0.001);
                    [ret,gstatus]=AndorGetStatus;
                    CheckWarning(ret);
                end
                if(ishghandle(h))
                    key = get(fig,'CurrentCharacter');
                    c_info = getCursorInfo(dcm_obj);
                    if Andor.ReadMode == 4 %image mode
                        [ret, imageData] = GetMostRecentImage(XPixels * YPixels);
                        CheckWarning(ret);
                        if ret == atmcd.DRV_SUCCESS
                            %display the acquired image
                            I=flipud(transpose(reshape(imageData, XPixels, YPixels)));
                            figure(1000);
                            set(h,'CData',I);axis image;
                            drawnow;
                        end
                    elseif Andor.ReadMode == 0 %FVB mode
                        [ret, img] = GetMostRecentImage(XPixels);
                        figure(1000);
                        h = plot(Andor.AxisWavelength, img.');axis([min(Andor.AxisWavelength) max(Andor.AxisWavelength) min(img(:)) max(img(:))]);
                        xlabel('Raman shifts (nm)'); ylabel('Raman counts');grid on;
                        CheckWarning(ret);
                        if strcmp(get(dcm_obj,'Enable'),'on')
                            if ~isempty(c_info) || ~isempty(old_info)
                                datacursors = dcm_obj.createDatatip(h);
                                if ~isempty(c_info)
                                    datacursors.Cursor.Position = [Andor.AxisWavelength(c_info.DataIndex), img(c_info.DataIndex) 0];
                                    old_info = c_info;
                                    c_info = [];
                                elseif ~isempty(old_info)
                                    datacursors.Cursor.Position = [Andor.AxisWavelength(old_info.DataIndex), img(old_info.DataIndex) 0];
                                end
                            end
                        end
                    end
                else
                    break;
                end
            end
            fprintf('Done!\n');
            if strcmp(key,'0')
                close(fig);
            end
            [ret]=SetShutter(1, 2, 1, 1); %close shutter
            CheckWarning(ret);
        end
        function imageData = acquire0(Andor) %acquire in scilence.
            [ret] = PrepareAcquisition();
            [ret]=SetShutter(1, 0, 0, 0); %set shutter open;
            CheckWarning(ret);
            [ret] = StartAcquisition();
            CheckWarning(ret);
            
            gstatus = 0;
            while(gstatus ~= atmcd.DRV_IDLE)
                [ret,gstatus]=AndorGetStatus;
                CheckWarning(ret);
            end
            acmod = Andor.AcquisitionMode;
            rdmod = Andor.ReadMode;
            ccdsize = Andor.CCD_Pixels;
            if   acmod == 1 || acmod == 2   %single scan and Accumulate
                if rdmod == 0 %FVB mode
                    [ret, imageData] = GetMostRecentImage(ccdsize(1));
                    CheckWarning(ret);
                elseif rdmod == 4 %Image mode;
                    [ret, imageData] = GetMostRecentImage(ccdsize(1) * ccdsize(2));
                    CheckWarning(ret);
                    if ret == atmcd.DRV_SUCCESS
                        imageData = flipud(transpose(reshape(imageData, ccdsize(1), ccdsize(2))));
                    end
                end
            elseif acmod == 3 %kinetics
                frameCount = Andor.NumberKinetics;
                if rdmod == 0 %FVB mode
                    imageData = zeros(frameCount,ccdsize(1));
                    for currentSeries = 1:frameCount
                        [ret, imageData(currentSeries,:)] = GetOldestImage(ccdsize(1));
                        CheckWarning(ret);
                    end
                elseif rdmod == 4 %Image mode;
                    imageData = zeros(frameCount,ccdsize(2),ccdsize(1));
                    for currentSeries = 1:frameCount
                        [ret, img] = GetMostRecentImage(ccdsize(1) * ccdsize(2));
                        CheckWarning(ret);
                        imageData(currentSeries,:,:) = flipud(transpose(reshape(img, ccdsize(1), ccdsize(2))));
                    end;
                end
            else
                fprintf('Acquisition mode %d is not supported now.',Andor.AcquisitionMode);
            end
        end
        
        function imageData = acquire(Andor)
            [ret,gstatus]=AndorGetStatus();
            CheckWarning(ret);
            if(gstatus ~= atmcd.DRV_IDLE)
                disp('System is busy, please try it again later.');
                return;
            end
            
            [ret]=SetShutter(1, 0, 0, 0); %set shutter open;
            CheckWarning(ret);
            
            [ret] = PrepareAcquisition();
            [ret] = StartAcquisition();
            CheckWarning(ret);
            
            fprintf('Acquiring.......');
            gstatus = 0;
            while(gstatus ~= atmcd.DRV_IDLE)
                [ret,gstatus]=AndorGetStatus;
                CheckWarning(ret);
            end
            fprintf('...done\n');
            
            acmod = Andor.AcquisitionMode;
            rdmod = Andor.ReadMode;
            ccdsize = Andor.CCD_Pixels;
            if   acmod == 1 || acmod == 2   %single scan and Accumulate
                if rdmod == 0 %FVB mode
                    [ret, imageData] = GetMostRecentImage(ccdsize(1));
                    CheckWarning(ret);
                elseif rdmod == 4 %Image mode;
                    [ret, imageData] = GetMostRecentImage(ccdsize(1) * ccdsize(2));
                    CheckWarning(ret);
                    if ret == atmcd.DRV_SUCCESS
                        imageData = flipud(transpose(reshape(imageData, ccdsize(1), ccdsize(2))));
                    end
                end
            elseif acmod == 3 %kinetics
                frameCount = Andor.NumberKinetics;
                if rdmod == 0 %FVB mode
                    imageData = zeros(frameCount,ccdsize(1));
                    for currentSeries = 1:frameCount
                        [ret, imageData(currentSeries,:)] = GetOldestImage(ccdsize(1));
                        CheckWarning(ret);
                    end
                elseif rdmod == 4 %Image mode;
                    imageData = zeros(frameCount,ccdsize(2),ccdsize(1));
                    for currentSeries = 1:frameCount
                        [ret, img] = GetMostRecentImage(ccdsize(1) * ccdsize(2));
                        CheckWarning(ret);
                        imageData(currentSeries,:,:) = flipud(transpose(reshape(img, ccdsize(1), ccdsize(2))));
                    end;
                end
            else
                fprintf('Acquisition mode %d is not supported now.',Andor.AcquisitionMode);
            end
            
            [ret]=SetShutter(1, 2, 1, 1); %set shutter automatic;
            CheckWarning(ret);
        end
        function imageData = getImage(Andor) %get images or spectra just acquired
            imageData = [];
            [ret,gstatus]=AndorGetStatus;
            CheckWarning(ret);
            if gstatus ~= atmcd.DRV_IDLE
                return;
            end
            acmod = Andor.AcquisitionMode;
            rdmod = Andor.ReadMode;
            [ret,XPixels, YPixels]=GetDetector;
                        CheckWarning(ret);
            if   acmod == 1 || acmod == 2   %single scan and Accumulate
                if rdmod == 4 %Image mode;
                    [ret] = StartAcquisition();
                    CheckWarning(ret);
                    
                    [ret] = WaitForAcquisition();
                    CheckWarning(ret);
                    [ret, imageData] = GetMostRecentImage(XPixels * YPixels);
                    CheckWarning(ret);
                    if ret == atmcd.DRV_SUCCESS
                        imageData = flipud(transpose(reshape(imageData, XPixels, YPixels)));
                    end
                    %                     imageData = flipud(transpose(reshape(imageData, ccdsize(1), ccdsize(2))));
                    
                end
            elseif acmod == 3 %kinetics
                frameCount = Andor.NumberKinetics;
                if rdmod == 4 %Image mode;
                    imageData = zeros(frameCount,YPixels,XPixels);
                    for currentSeries = 1:frameCount
                        [ret, img] = GetMostRecentImage(XPixels * YPixels);
                        CheckWarning(ret);
                        imageData(currentSeries,:,:) = flipud(transpose(reshape(img, XPixels, YPixels)));
                    end;
                end
            else
                fprintf('Acquisition mode %d is not supported now.',Andor.AcquisitionMode);
            end
        end
    end
    
    methods (Access = protected)
        
    end
    
    methods (Access = private)
        function CropHeight = setCropParam(Andor,newCropHeight)
            currentMode = Andor.ReadMode;
            if currentMode~=0; % only under FVB mode can set crop;
                Andor.ReadMode = 0; %Set to FVB
            end
            if newCropHeight<=0 || newCropHeight >= Andor.CCDsize(2)
                newCropHeight = Andor.CCDsize(2);
                [ret] = SetCropMode(0, newCropHeight, 0);
            else
                [ret] = SetCropMode(1, newCropHeight, 0);
            end
            CropHeight = newCropHeight;
            CheckWarning(ret);
            if currentMode~=0; % only under FVB mode can set crop;
                Andor.ReadMode = currentMode; %reset read mode;
            end
        end
        function noGrating = GetNumberofAvailableGratings(Andor)
            noGrating = 0;
            [ret, present] = ShamrockGratingIsPresent(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if present
                [ret, noGrating] = ShamrockGetNumberGratings(Andor.ShamrockDev);
                ShamrockCheckWarning(ret);
            end
        end
        function AxisWavelength = GetAxisWavelength(Andor)
            [ret, NumberPixels] = ShamrockGetNumberPixels(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            [ret, AxisWavelength] = ShamrockGetCalibration(Andor.ShamrockDev, NumberPixels);
            ShamrockCheckWarning(ret);
        end
        function CentralWavelength = GetCentralWavelength(Andor)
            [ret, present] = ShamrockWavelengthIsPresent(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if present
                [ret, CentralWavelength] = ShamrockGetWavelength(Andor.ShamrockDev);
                ShamrockCheckWarning(ret);
            end
        end
        function CentralWavelength = SetCentralWavelength(Andor,newCentralWavelength)
            [ret, present] = ShamrockWavelengthIsPresent(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if present
                [ret] = ShamrockSetWavelength(Andor.ShamrockDev, newCentralWavelength);
                ShamrockCheckWarning(ret);
                CentralWavelength = newCentralWavelength;
                Andor.AxisWavelength = GetAxisWavelength(Andor);
                fprintf('Central wavelength has been set to %fnm\n',Andor.CentralWavelength);
            end
        end
        function SlitWidth = GetSlitWidth(Andor)
            [ret, present] = ShamrockSlitIsPresent(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if present
                [ret, SlitWidth] = ShamrockGetSlit(Andor.ShamrockDev);
                ShamrockCheckWarning(ret);
            end
        end
        function SlitWidth = SetSlitWidth(Andor,newSlitWidth)
            [ret, present] = ShamrockSlitIsPresent(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if present
                [ret] = ShamrockSetSlit(Andor.ShamrockDev,newSlitWidth);
                ShamrockCheckWarning(ret);
                if ret == 20202 %successful;
                    SlitWidth = newSlitWidth;
                end
            end
        end
        function CurrentGrating = GetCurrentGrating(Andor)
            [ret, grating] = ShamrockGetGrating(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if ret == 20202 %successful;
                CurrentGrating = grating;
            end
        end
        function CurrentGrating = SetCurrentGrating(Andor,newGrating)
            [ret, grating] = ShamrockGetGrating(Andor.ShamrockDev);
            ShamrockCheckWarning(ret);
            if newGrating > 0 && newGrating <= Andor.noGratings && grating~=newGrating
                [ret] = ShamrockSetGrating(Andor.ShamrockDev, newGrating);
                if ret == 20202
                    CurrentGrating = newGrating;
                end
                [ret, Lines, Blaze, Home, Offset] = ShamrockGetGratingInfo(Andor.ShamrockDev, Andor.CurrentGrating);
                ShamrockCheckWarning(ret);
                fprintf('New grating Info: %d lines/mm, Blaze@%snm, Offset: %d, SlitWidth: %dum\n',Lines,Blaze,Offset,Andor.SlitWidth);
            end
        end
        
    end
    
    methods (Static = true)
        function status = isAndorIdle()
            [ret,gstatus]=AndorGetStatus;
            CheckWarning(ret);
            if gstatus == atmcd.DRV_IDLE
                status = 1;
            else
                status = 0;
            end
        end
        function status = startAcquire() % send a command to start the acqistion;
            [ret,gstatus]=AndorGetStatus();
            CheckWarning(ret);
            if(gstatus ~= atmcd.DRV_IDLE)
                status = 0;
                return;
            end
            status = 1;
            [ret]=SetShutter(1, 0, 0, 0); %set shutter open;
            CheckWarning(ret);
            [ret] = PrepareAcquisition();
            CheckWarning(ret);
            [ret] = StartAcquisition();
            CheckWarning(ret);
        end
        function abortAcquire()
            ret = AbortAcquisition();
            CheckWarning(ret);
        end
        function releaseAndor()
            [ret]=SetShutter(1, 2, 1, 1); %close shutter
            CheckWarning(ret);
            [ret, iCoolerStatus] = IsCoolerOn();
            CheckWarning(ret);
            if iCoolerStatus
                ret = CoolerOFF();
                CheckWarning(ret);
            end
            [ret] = SetCoolerMode(1); % keep cooler on;
            CheckWarning(ret);
            
            [ret]=AndorShutDown;
            CheckWarning(ret);
            [ret] = ShamrockClose();
            ShamrockCheckWarning(ret);
        end
    end
end

function txt = myupdatefcn(~,event_obj)
pos = get(event_obj,'Position');
Ind = get(event_obj, 'DataIndex');
if ~isempty(pos) && ~isempty(Ind)
    txt = {sprintf('Wavelength: %5.2fnm',pos(1)),sprintf('Raman Intensity: %5.0f',pos(2)),sprintf('Data Index: %d',Ind(1))};
end
end