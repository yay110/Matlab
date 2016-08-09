classdef finesse < handle
      
    %This class is for control the finesse laser parameters, such as power
    %and shutter.
    
    % script is written by Zhengyi Yang (zy6@st-andrews.ac.uk)
    % Created on 09/08/2016
    
    properties
        port;
        power;
        shutter;
        target;
    end
    
    methods
        function laser = finesse(port)
            if nargin < 1
                laser.port = 'COM1';
            else
                laser.port = port;
            end
            laser.target = serial(laser.port);
            set(laser.target,'BaudRate',19200);
            laser.target.Terminator = 'CR';
            fopen(laser.target);
        end
        
        function laser = open(laser)
            fprintf(laser.target,'SHUTTER OPEN');
            laser.shutter = 'Open';
        end
        
        function laser = close(laser)
            fprintf(laser.target, 'Shutter Close');
            laser.shutter = 'Closed';
        end
        
        function laser = set.power(laser,n)
            fprintf(laser.target,strcat('POWER=',num2str(n)));
%             fprintf(laser.target,'POWER?');
%             response = fscanf(laser.target);
%             laser.power = str2double(response(1:6));
        end
        
        function power = get.power(laser)
            fprintf(laser.target,'POWER?');
            response = fscanf(laser.target);
            power = str2double(response(1:6));
        end
        
        function laser = clear(laser)
            fclose(laser.target);
            delete(laser.target);
        end
    end
end