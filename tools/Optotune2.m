classdef Optotune2 < handle
    
    % Optotune class is designed to control electro tunable lens from optotune.com
    % The script is written by R. Spesyvtsev from the Optical manipulation group (http://www.st-andrews.ac.uk/~photon/manipulation/)
    % Last modification by Zhengyi Yang
    %
    % constructor Optotune(com port), example lens = Optotune('COM9');
    % connect: lens = lens.Open();
    % set current to 100 mA:  lens = lens.setCurrent(100);
    % lens can be set to different modes:
    % Current/Focal Power/Analog/Sinusoidal/Rectangular/Triangular
    % mode LowerCurrent/UpperCurrent/Frequency includes the parameters.
    
    % checkError will check the response if there is error and what error
    
    properties
        etl_port;
        port;
        status;
        response;
        
        temperature = 0;
        current = 0; %% in mAmpers
        max_current = 0; % in mAmpers
        min_current = 0; % in mAmpers
        calibration = 1;  % in mAmpers/micrometer
        
        modeLowerCurrent = 0;
        modeUpperCurrent = 0;
        modeFrequency = 1;
        
        max_bin = 0;
        time_pause = 0.3;
        time_laps = 0.01;
        last_time_laps;
        
        error = false;
    end
    
    properties (Dependent)
        temperatureLimits;
        mode;
        analogInput;
    end
    
    methods
        function lens=Optotune2(port)
            if (nargin<1)
                lens.port='COM3';
            else
                lens.port = port;
            end
        end
        
        function lens = Open(lens)
            %% Setting up initial parameters for the com port
            lens.etl_port = serial(lens.port);
            lens.etl_port.Baudrate=115200;
            lens.etl_port.StopBits=1;
            lens.etl_port.Parity='none';
            
            fopen(lens.etl_port);
            lens.status = lens.etl_port.Status;   %%% checking if initialization was completed the status should be "open"
            
            %% Initialize communication
            fprintf(lens.etl_port, 'Start'); %%% initiating the communication
            lens.last_time_laps = checkStatus(lens);
            %lens.etl_port.BytesAvailable;  %%% checking number of bytes in the response
            fscanf(lens.etl_port);  %%% reading out the response which should read "Ready"
            if lens.etl_port.BytesAvailable
                fread(lens.etl_port,lens.etl_port.BytesAvailable);
            end
            
            %get initial information about UpperCurrentLimte and MaxBin;
            lens.max_bin();
            lens.max_current();
            %Get initial parameters
            lens.modeFrequency;
            lens.modeLowerCurrent;
            lens.modeUpperCurrent;
            lens.temperature;
            
        end
        
        function temperature = get.temperature(lens)
            command = append_crc('TA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            %pause(lens.time_pause);
            %lens.etl_port.BytesAvailable;  %%% checking number of bytes in the response
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            x = lens.response(4)*(hex2dec('ff')+1) + lens.response(5);
            temperature = x*0.0625;
        end
        
        function current = get.current(lens)
            %% Get current in mA %%%%%%%%
            command = append_crc(['Ar'-0 0 0]);
            fwrite(lens.etl_port, command);
            %pause(lens.time_pause);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            current = lens.response(2)*(hex2dec('ff')+1) + lens.response(3);
            current = current * lens.max_current / (lens.max_bin+1);
        end
        
        function lens = set.current(lens, ci) %% Set current in mA via ci variable
            %% Set current %%%%%%%%
            set_i = (floor(ci*(lens.max_bin+1) / lens.max_current));
            LB = mod(set_i,256); %% low byte
            HB = (set_i-LB)/256; %% high byte
            command = append_crc(['Aw'-0 HB LB]);
            fwrite(lens.etl_port, command);
            %pause(lens.time_pause);
        end
        
        %% setting temperature limits for operation in focal power mode
        function lens = set.temperatureLimits(lens,~)
            % The function has not been finished
            temp_high = 80/0.0625;  %% Upper limit 80 degrees (see Optotune manual for available limits)
            temp_low = 0/0.0625;     %% Lower limit 0 degrees
            LBH = mod(temp_high,256); %% low byte
            HBH = (temp_high-LBH)/256; %% high byte
            LBL = mod(temp_low,256); %% low byte
            HBL = (temp_low-LBL)/256; %% high byte
            command = append_crc(['PrTA'-0 HBH LBH HBL LBL]);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
        end
        
        function max_bin = get.max_bin(lens)
            command = append_crc(['CrUA'-0 0 0]);
            fwrite(lens.etl_port, command);
            %pause(lens.time_pause);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            max_bin = lens.response(4)*(hex2dec('ff')+1) + lens.response(5)+1;  %% hardware current limit usually 4095;
        end
        
        %% Get current limit %%%%%%%%
        function max_current = get.max_current(lens)
            command = append_crc(['CrMA'-0 0 0]);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            %pause(lens.time_pause);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            max_current = lens.response(4)*(hex2dec('ff')+1) + lens.response(5);  %% software current limit usually 292.84 mA;
            max_current = max_current / 100; %% reads current in mili ampers
        end
        
        %% set the lens to different mode
        function lens = currentMode(lens)
            command = append_crc('MwDA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MDA')
                logMessage('Lens set to Current Mode succesfully');
            else
                checkError(lens);
            end
        end
        
        function lens = focalPowerMode(lens)
            command = append_crc('MwCA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MCA')
                logMessage('Lens set to Focal Power Mode succesfully');
            else
                checkError(lens);
            end
        end
        
        function lens = analogMode(lens)
            command = append_crc('MwAA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MAA')
                logMessage('Lens set to Analog Mode succesfully');
            else
                checkError(lens);
            end
        end
        
        function lens = sinusoidalMode(lens)
            command = append_crc('MwSA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MSA')
                logMessage('Lens set to Sinusoidal Signal Mode succesfully');
            else
                checkError(lens);
            end
%             lens = lens.setModeFrequency(lens.modeFrequency);
%             lens = lens.setModeLowerCurrent(lens.modeLowerCurrent);
%             lens = lens.setModeUpperCurrent(lens.modeUpperCurrent);
        end
        
        function lens = rectangularMode(lens)
            command = append_crc('MwQA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MQA')
                logMessage('Lens set to Rectangular Signal Mode succesfully');
            else
                checkError(lens);
            end
%             lens = lens.setModeFrequency(lens.modeFrequency);
%             lens = lens.setModeLowerCurrent(lens.modeLowerCurrent);
%             lens = lens.setModeUpperCurrent(lens.modeUpperCurrent);
        end
        
        function lens = triangularMode(lens)
            command = append_crc('MwTA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:3))'), 'MTA')
                logMessage('Lens set to Triangular Signal Mode succesfully');
            else
                checkError(lens);
            end
%             lens = lens.setModeFrequency(lens.modeFrequency);
%             lens = lens.setModeLowerCurrent(lens.modeLowerCurrent);
%             lens = lens.setModeUpperCurrent(lens.modeUpperCurrent);
        end
        
        %% get the current mode the lens is running at
        function mode = get.mode(lens)
            command = append_crc('MMA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            checkError(lens);
            switch lens.response(4)
                case 1
                    logMessage('Lens is driven in Current mode');
                case 2
                    logMessage('Lens is driven in Sinusoidal Singal mode');
                case 3
                    logMessage('Lens is driven in Triangular mode');
                case 4
                    logMessage('Lens is driven in Retangular mode');
                case 5
                    logMessage('Lens is driven in FocalPower mode');
                case 6
                    logMessage('Lens is driven in Analog mode');
                case 7
                    logMessage('Lens is driven in Position Controlled Mode');
            end
        end
        
        %% set parameters for the mode controlling
        
        % set signal generator upper current limit
        function lens = set.modeUpperCurrent(lens,ci)
            if ci > lens.max_current || ci < 0
                logMessage('The current should be between 0 and %.2f mA', lens.max_current);
                ci = lens.max_current;
            else
                if ci < lens.modeLowerCurrent
                    lens = lens.setModeLowerCurrent(ci);
                end
            end
            set_i = (floor(ci*(lens.max_bin+1) / lens.max_current));
            LB = mod(set_i,256); %% low byte
            HB = (set_i-LB)/256; %% high byte
            command = append_crc(['PwUA'-0 HB LB 0 0]);
            fwrite(lens.etl_port, command);
            checkError(lens);
            if lens.error == true
                logMessage('Error occurred when setting modeUpperCurrent');
            end
            lens.modeUpperCurrent;
        end
        
        % get signal generator upper current limit
        function modeUpperCurrent = get.modeUpperCurrent(lens)
            command = append_crc(['PrUA'-0 0 0 0 0]);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:4))'), 'MPAU')
                modeUpperCurrent = (lens.response(5)*(hex2dec('ff')+1) + lens.response(6))*lens.max_current/lens.max_bin;
            else
                modeUpperCurrent = 0;
                checkError(lens);
                if lens.error == true
                    logMessage('Error occurred when getting modeUpperCurrent');
                end
            end
        end
        
        % set signal generator lower current limit
        function lens = set.modeLowerCurrent(lens,ci)
            if ci > lens.max_current || ci < 0
                logMessage('The current should be between 0 and %.2f mA', lens.max_current);
                ci = 0;
            else
                if ci > lens.modeUpperCurrent
                    lens = lens.setModeUpperCurrent(ci);
                end
            end
            set_i = (floor(ci*(lens.max_bin+1) / lens.max_current));
            LB = mod(set_i,256); %% low byte
            HB = (set_i-LB)/256; %% high byte
            command = append_crc(['PwLA'-0 HB LB 0 0]);
            fwrite(lens.etl_port, command);
            checkError(lens);
            if lens.error == true
                logMessage('Error occurred when setting modeLowerCurrent');
            end
            lens.modeLowerCurrent;
        end
        
        % get signal generator lower current limit
        function modeLowerCurrent = get.modeLowerCurrent(lens)
            command = append_crc(['PrLA'-0 0 0 0 0]);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:4))'), 'MPAL')
                modeLowerCurrent = (lens.response(5)*(hex2dec('ff')+1) + lens.response(6))*lens.max_current/lens.max_bin;  %% software current limit usually 292.84 mA;
            else
                modeLowerCurrent = 0;
                checkError(lens);
                if lens.error == true
                    logMessage('Error occurred when getting modeLowerCurrent');
                end
            end
        end
        
        % set signal generator lower current limit
        function lens = set.modeFrequency(lens,ci)
            if ci < 0.1
                logMessage('The minimum frequency is 0.1 Hz!');
                ci = 0.1;
            end
            set_i = ci*1000;
            B4 = mod(set_i,256); %% 4th byte
            B3 = mod((set_i-B4)/256,256); %% third byte
            B2 = mod((set_i-B3*256-B4),256); %% second byte
            B1 = mod((set_i-B2*2^16-B3*256-B4),256); %% first byte
            command = append_crc(['PwFA'-0 B1 B2 B3 B4]);
            fwrite(lens.etl_port, command);
            checkError(lens);
            if lens.error == true
                logMessage('Error occurred when setting modeFrequency');
            end
            lens.modeFrequency;
        end
        
        % get signal generator lower current limit
        function modeFrequency = get.modeFrequency(lens)
            command = append_crc(['PrFA'-0 0 0 0 0]);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            if strcmp(cellstr(char(lens.response(1:4))'), 'MPAF')
                modeFrequency = (lens.response(5) * 2^24 + lens.response(6) * 2^16 + lens.response(7) * 2^8 + lens.response(8))/1000;  %% software current limit usually 292.84 mA;
            else
                checkError(lens);
                if lens.error == true
                    logMessage('Error occurred when getting modeFrequency');
                end
            end
        end
        
        
        function analogInput = get.analogInput(lens)
            command = append_crc('GAA'-0);
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            analogInput = lens.response(4)*(hex2dec('ff')+1) + lens.response(5);
        end
        
        function lens = testError(lens)
            command = 'St4rt';
            fwrite(lens.etl_port, command);
            lens.last_time_laps = checkStatus(lens);
            lens.response = fread(lens.etl_port,lens.etl_port.BytesAvailable);
            checkError(lens);
        end
        
        %%  Closing the port when finished using it %%%%%%%%%%%%
        function lens = Close(lens)
            lens = lens.currentMode();
            lens.current = 0;
            fclose(lens.etl_port);
            delete(lens.etl_port);
            clear lens.etl_port
            
            lens.status = 'closed';
            lens.response = 'Shut down';
            
        end
        
        function tElapsed = checkStatus(lens)
            bts = lens.etl_port.BytesAvailable;  %%% checking number of bytes in the response
            tStart = tic;
            tElapsed = 0;
            while (bts ==0) || (tElapsed >5)
                bts = lens.etl_port.BytesAvailable;  %%% checking number of bytes in the response
                pause(lens.time_laps);
                tElapsed = toc(tStart);
            end
        end
        
        %Check if there is error message and identify what it is, display.
        function checkError(lens)
            lens.error = false;
            if char(lens.response(1)) == 'E'
                switch lens.response(2)
                    case 1
                        logMessage('CRC failed.');
                    case 2
                        logMessage('Command not available in firmware type.');
                    case 3
                        logMessage('Command not regongnized.');
                    case 4
                        logMessage('Lens is not compatible with firmware.');
                        
                end
                lens.error = true;
            end
        end
    end
end
