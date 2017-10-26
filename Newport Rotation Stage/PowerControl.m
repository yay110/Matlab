cd 'C:\Users\OMG\Documents\GitHub\Matlab\Newport Rotation Stage';
NewportESP = serial('COM5');
set(NewportESP,'BaudRate',19200);
NewportESP.Parity = 'none';
NewportESP.Terminator = 'CR/LF';
smcAxis =1;

wavelength = 800;

file = dir(['*doubleWP-PowerCalibration' num2str(wavelength) 'nm*.mat']);
load(file.name);
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);

fclose(NewportESP);
fopen(NewportESP);



%deg = rot_calibr(0.5,f1);smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(deg));fprintf(laser,smcErrCmnd); % set absolute position

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(min(power(:)),f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.05,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.1,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.2,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.3,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.5,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(1.0,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(1.5,f1))));

fprintf(NewportESP,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(max(power(:)),f1))));

fprintf(NewportESP,'1PA0');

fclose(NewportESP);
delete(NewportESP);
laser.Close;
