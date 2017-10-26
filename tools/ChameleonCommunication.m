%% Instrument Connection

% Find a serial port object.
Chameleon = instrfind('Type', 'serial', 'Port', 'COM1', 'Tag', '');

% Create the serial port object if it does not exist
% otherwise use the object that was found.
if isempty(Chameleon)
    Chameleon = serial('COM1');
else
    fclose(Chameleon);
    Chameleon = Chameleon(1);
end

% Connect to instrument object, obj1.
fopen(Chameleon);

set(Chameleon, 'BaudRate', 19200);
set(Chameleon, 'Terminator', {'CR/LF','CR/LF'});


fprintf(Chameleon,'LASER=1');
meas = fscanf(Chameleon)


fprintf(Chameleon,'?LIGHT');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)

fprintf(Chameleon,'SHUTTER=0');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)
fprintf(Chameleon,'?S');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)

fprintf(Chameleon,'?VW');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)

%wavelength Tunning
fprintf(Chameleon,'WAVELENGTH=795');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)
fprintf(Chameleon,'?TS');
Chameleon.BytesAvailable
meas = fscanf(Chameleon)


fclose(Chameleon);
