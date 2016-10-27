classdef GilsonPump < handle
    
    %This class is for control the Gilson Minipuls 3 pump remotely, through
    %a Arduino board. At back of the pump, there
    
    % script is written by Zhengyi Yang (zy6@st-andrews.ac.uk)
    % Created on 09/08/2016
    
    properties
        port;
        target;
    end
    
    methods
        function pump = GilsonPump(port)
            if nargin < 1
                pump.port = 'COM17';
            else
                pump.port = port;
            end
            pump.target = serial(pump.port);
        end
        
        function pump = open(pump)
            fopen(pump.target);
            pause(2);
            fprintf(pump.target,'A');
        end
        
        function pump = on(pump)
            fprintf(pump.target,'B');
        end
        
        function pump = off(pump)
            fprintf(pump.target,'A');
        end
        
        function pump = clear(pump)
            fclose(pump.target);
        end
    end
    
end