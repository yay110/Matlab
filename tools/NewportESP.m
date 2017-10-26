classdef NewportESP < handle
    
    properties
        SerialPort;
        Axis;
        absPosition;
    end
    
    methods
        function NewportESP = NewportESP(Port,Axis)
            if nargin<1
                Port = 'COM5';
            end
            
            if nargin<2
                Axis = 1;
            end
            NewportESP.Axis = Axis;
            
            SerialPort = serial(Port);
            NewportESP.SerialPort = SerialPort;
            % Connect to instrument object, obj1.
            set(SerialPort, 'BaudRate', 19200);
            set(SerialPort, 'Parity', 'none');
            set(SerialPort, 'Terminator', {'CR/LF','CR/LF'});
            fclose(SerialPort);
            fopen(SerialPort);
            disp('Connected to NewportESP......');
        end
        
        function absPosition = get.absPosition(NewportESP)
            axis = NewportESP.Axis;
            Cmnd = [num2str(axis),'TP'];
            absPosition = str2double(query(NewportESP.SerialPort,Cmnd)); %% check absolute position
        end
        
        function NewportESP = set.absPosition(NewportESP,absPosition)
            Cmnd = strcat(num2str(NewportESP.Axis),'PA',num2str(absPosition));
            fprintf(NewportESP.SerialPort,Cmnd); % set absolute position
        end
        
        function NewportESP = Close(NewportESP)
           fclose(NewportESP.SerialPort);
        end
    end
    
    
end