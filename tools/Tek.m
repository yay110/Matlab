classdef Tek < handle
    
    % FunGenTek class is designed to control the Tektronix AFG1000 seris
    % function generator. AFG1022 is the model we have.
    
    % The script is written by Zhengyi Yang from Optical manipulation group
    % St Andrews, UK
    % Last edited on 14/07/2016
    
    properties
        target;
        output;
        mode;
        frequency;
        amplitude;      %Unit V;
        offset;         %Unit V;
        phaseOffset;    %Unit deg
        burstCycles;
        burst;
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
            
            Tek.mode = 1;
        end
        
        % basic function such as choosing mode, change frequency, amplitude
        % and the offset
        function Tek = set.mode(Tek,n)
            switch n
                case 1
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE SIN');
                    Tek.mode = 'Sin';
                case 2
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE SQUARE');
                    Tek.mode = 'Squre';
                case 3
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE Ramp');
                    Tek.mode = 'Ramp';
                case 4
                    fprintf(Tek.target, 'Source1:FUNCTION:SHAPE Pulse');
                    Tek.mode = 'Pulse';
            end
        end
        
        function Tek = set.frequency(Tek,n)
            fprintf(Tek.target, ['SOURCE1:FREQUENCY ', num2str(n), 'HZ']);
            Tek.frequency = n;
        end
        
        function Tek = set.amplitude(Tek,n)
            fprintf(Tek.target,['SOURCE1:VOLTAGE:AMPLITUDE ',num2str(n),'VPP']);
            Tek.amplitude = n;
        end
        
        function Tek = set.offset(Tek,n)
            fprintf(Tek.target,['SOURCE1:VOLTAGE:OFFSET ',num2str(n),'V']);
            Tek.offset = n;
        end
        
        % Output status of source1
        function Tek = outputOn(Tek)
            fprintf(Tek.target, 'OUTPUT1:STATE ON');
            Tek.output = 'On';
        end
        
        function Tek = outputOff(Tek)
            fprintf(Tek.target, 'OUTPUT1:STATE OFF');
            Tek.output = 'Off';
        end
        
        %Burst mode related function
        function Tek = burstOn(Tek)
            fprintf(Tek.target, 'SOURCE1:BURST:STATE ON');
                      Tek.burst = 'On';
        end
        
        function Tek = burstOff(Tek)
            fprintf(Tek.target, 'SOURCE1:BURST:STATE OFF');
                      Tek.burst = 'Off';
        end
        
        function Tek = set.burstCycles(Tek,n)
            fprintf(Tek.target, ['SOURCE1:BURST:ncycles ',num2str(n)]);
            Tek.burstCycles = n;
        end
        
        %important function needed!!!!!
        % change the relative phase for the burst mode
        % this function is not added because there is no mention of it on
        % the programmer manual.
        % need to contact the company to add it. 
        
        
        %% This function is needed to scan the function generator according
        %  to tunable lens
        % The values needs to be calibrated for each system, for each
        % frequency.
        function Tek = scan(Tek,n)
            Tek.burstOn;
            Tek.frequency = n;
            Tek.burstCycles = 1;
            switch n
                case 10
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                    Tek.burstCycles = 50;
                case 20
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 25
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 30
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 50
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 80
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 200
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 250
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 300
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                case 320
                    Tek.offset = 0.106;                    
                    Tek.amplitude = 0.073;
                otherwise
                    disp('This frequency is not calibrated!');
            end
            Tek.outputOn;
        end
        
    end
    
end
