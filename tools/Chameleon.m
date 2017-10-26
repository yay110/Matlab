classdef Chameleon < handle
    
    %Chameleon Class for remote control of Chameleon Ultra II from Coherent
    properties
        SerialPort;
        Wavelength;
        WavelengthOffSet = -5;
        TuningStatus;
        StatusShutter = 0;
    end
    
    methods
        function Chameleon = Chameleon(Port)
            if nargin<1
                Port = 'COM1';
            end
            
            SerialPort = serial(Port);
            Chameleon.SerialPort = SerialPort;
            fopen(SerialPort);
            set(SerialPort, 'BaudRate', 19200);
            set(SerialPort, 'Terminator', {'CR/LF','CR/LF'});
        end
        
        function Chameleon = SOpen(Chameleon)
            disp(query(Chameleon.SerialPort,'SHUTTER=1'));
            Chameleon.StatusShutter = 1;
        end
        function Chameleon = SClose(Chameleon)
            disp(query(Chameleon.SerialPort,'SHUTTER=0'));
            Chameleon.StatusShutter = 0;
        end
        
        function Chameleon = set.Wavelength(Chameleon,Wavelength)
            WavelengthOffset = Chameleon.WavelengthOffSet;
            fprintf(Chameleon.SerialPort,strcat('WAVELENGTH=',num2str(Wavelength+WavelengthOffset)));
            delay(5);
            [~] = fscanf(Chameleon.SerialPort);
            meas = query(Chameleon.SerialPort,'?VW');
            Wavelength = str2double(meas(16:end))-WavelengthOffset;
            disp(['Now the wavelength is tuned to ' num2str(Wavelength) 'nm.']);
        end
        
        function Wavelength = get.Wavelength(Chameleon)
            WavelengthOffset = Chameleon.WavelengthOffSet;
            meas = query(Chameleon.SerialPort,'?VW');
            Wavelength = str2double(meas(16:end))-WavelengthOffset;
        end
        
        function Chameleon = Close(Chameleon)
            Chameleon.SClose;
            fclose(Chameleon.SerialPort);
            disp('Chameleon is disconnected.')
        end
        
    end
end