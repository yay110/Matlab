% This program shows the basic command to connect and operate the function
% generator Keithley 3390. 

% The 3390 unit is USBTMC (USB Test & Measurement Class Device). This is a
% Protocol which is built on top of USB that allows GPIB-like communication
% with USB devices. Thus you can talk to the 3390 as if it was connected
% via GPIB using the full SCPI command set as defined in the 3390 user
% manual. 
% To talk to the USBTMC instrument, you will first need to install
% NI-VISA(download and install the latest Keithley I/O layer, you will get
% NI-VISA runtime). Connect the device to computer via USB, then use MATLAB
% Test & Measurement Tool to test the equipment. 

% Find a VISA-USB object.
obj1 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x05E6::0x3390::1242585::0::INSTR', 'Tag', '');

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('NI', 'USB0::0x05E6::0x3390::1242585::0::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Connect to instrument object, obj1.
fopen(obj1);

%% Communicating with instrument object, obj1.

%Sin function with the amplitude and offset
fprintf(obj1, 'APPL:SIN 5 KHZ, 3.0 VPP, -2.5 V');
% DC with 0.3V amplitude
fprintf(obj1, 'APPL:DC 10 KHZ,  3 VPP, 0.3 V');
% Ramp function with amplitude 3V, offset -2.5V and 100% symmetry
fprintf(obj1, 'APPL:RAMP 5 KHZ, 3.0 VPP, -2.5 V');
fprintf(obj1, 'FUNCTION:RAMP:SYMMETRY 100');

fprintf(obj1, 'OUTPUT OFF');



% Disconnect from instrument object, obj1.
fclose(obj1);