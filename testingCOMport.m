s= serial('COM1');
set(s,'BaudRate',19200);
s.Terminator = 'CR';
fopen(s);
fprintf(s,'SHUTTER Close');
fprintf(s,'POWER?');
out = fscanf(s);
fprintf(s,'POWER=0.5');
out = fscanf(s);
fclose(s);
delete(s);
clear s