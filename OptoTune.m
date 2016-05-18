%Control the electrical  tuable lens from a vitual COM port
%Configuration setting are?
%Baudrate: 115200
%Parity: None
%Stop Bits: One
%Created by Zhengyi Yang(zy6@st-andrews.ac.uk) on 13th May 2016
%Last edited on 13th May 2016

%set the ETL to the right COM port
ETL=serial('COM3');
%set(ETL,'BaudRate',115200);
ETL.BaudRate = 115200;
fopen(ETL);

%Give the port a command
%Write and read data — Write data to the device using the fprintf or fwrite function, and read data from the device using the fgetl, fgets, fread, fscanf, or readasync function.
fprintf(ETL,'AW');
%get the reply from the COM port
reply = fscanf(ETL);

fclose(ETL);
delete(ETL);
clear ETL;