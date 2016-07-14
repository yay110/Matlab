classdef FunGenTek
    
    % FunGenTek class is designed to control the Tektronix AFG1000 seris
    % function generator. AFG1022 is the model we have.
    
    % The script is written by Zhengyi Yang from Optical manipulation group
    % St Andrews, UK
    % Last edited on 14/07/2016
    
    properties
        obj;
        mode = 'sin';
        frequency = 10;
        amplitude = 0;      %Unit V;
        offset    = 0;      %Unit V;
        phaseOffset = 90;   %Unit degree
        
        
    end
    
    methods
        function FunGenTek = FunGenTek()
            if (nargin <1)
                
                % Find a VISA-USB object.
                FunGenTek.obj = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0353::1618156::0::INSTR', 'Tag', '');
                
                % Create the VISA-USB object if it does not exist
                % otherwise use the object that was found.
                if isempty(FunGenTek.obj)
                    FunGenTek.obj = visa('NI', 'USB0::0x0699::0x0353::1618156::0::INSTR');
                else
                    fclose(FunGenTek.obj);
                    FunGenTek.obj = FunGenTek.obj(1);
                end
                
                
                % Connect to instrument object, obj1.
                fopen(FunGenTek.obj);
            end
        end
        
        function obj = setFrequency(obj,n)
            obj.obj.frequency = n;
            fprintf(obj.obj, ['SOURCE1:FREQUENCY ', n, 'HZ']);
        end
        
        
        
    end
    
end
