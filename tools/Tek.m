classdef Tek
    
    % FunGenTek class is designed to control the Tektronix AFG1000 seris
    % function generator. AFG1022 is the model we have.
    
    % The script is written by Zhengyi Yang from Optical manipulation group
    % St Andrews, UK
    % Last edited on 14/07/2016
    
    properties
        target;
        output;
        mode = 1;
        frequency = 10;
        amplitude = 0;      %Unit V;
        offset    = 0;      %Unit V;
        phaseOffset = 90;   %Unit degree
        
        
    end
    
    methods
        % initial the function generator and connect in Matlab
        function Tek = Tek()
            % Find a VISA-USB object.
            Tek.target = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0353::1618156::0::INSTR', 'Tag', '');
            % Create the VISA-USB object if it does not exist
            % otherwise use the object that was found.
            if isempty(Tek.target)
                Tek.target = visa('NI', 'USB0::0x0699::0x0353::1618156::0::INSTR');
            else
                fclose(Tek.target);
                Tek.target = Tek.target(1);
            end
            % Connect to instrument object, obj1.
            fopen(Tek.target);
            
        end
        
        % basic function such as choosing mode, change frequency, amplitude
        % and the offset
        function Tek = set.mode(Tek,n)
            switch n
                case 1
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE SIN');
                case 2
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE SQUARE');
                case 3
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE Ramp');
                case 4
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE Pulse');
            end
        end
        
        function Tek = set.frequency(Tek,n)
            fprintf(Tek.target, ['SOURCE1:FREQUENCY ', num2str(n), 'HZ']);
        end
        
        function Tek = set.amplitude(Tek,n)
            fprintf(Tek.target,['SOURCE1:VOLTAGE:AMPLITUDE ',num2str(n),'VPP']);
        end
        
        function Tek = set.offset(Tek,n)
            fprintf(Tek.target,['SOURCE1:VOLTAGE:OFFSET ',num2str(n),'V']);
        end
        
        % Output status of source1
        function Tek = outputOn(Tek)
            fprintf(Tek.target, 'OUTPUT1:STATE ON');
        end
        
        function Tek = outputOff(Tek)
            fprintf(Tek.target, 'OUTPUT1:STATE OFF');
        end
        
        %Burst mode related function
        
        
    end
    
end
