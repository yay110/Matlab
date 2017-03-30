classdef Arduino < handle
    
    %This class is for control the Gilson Minipuls 3 pump remotely, through
    %a Arduino board. At back of the pump, there
    
    % script is written by Zhengyi Yang (zy6@st-andrews.ac.uk)
    % Created on 09/08/2016
    
    properties
        port;
        target;
    end
    
    methods
        function board = Arduino(port)
            if nargin < 1
                board.port = 'COM17';
            else
                board.port = port;
            end
            board.target = serial(board.port);
        end
        
        function board = connect(board)
            fopen(board.target);
            pause(2);
            %pump off
            fprintf(board.target,'A');
            %laser off
            fprintf(board.target,'D');
        end
        
        function board = pump_on(board)
            fprintf(board.target,'B');
        end
        
        function board = pump_off(board)
            fprintf(board.target,'A');
        end
        
        function board = laser_on(board)
            fprintf(board.target,'C');
        end
        
        function board = laser_off(board)
            fprintf(board.target,'D');
        end
        
        function board = disconnect(board)
            fclose(board.target);
        end
    end
    
end