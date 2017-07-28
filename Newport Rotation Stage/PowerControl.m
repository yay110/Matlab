
esp_port = serial('COM11');
set(esp_port,'BaudRate',19200);
esp_port.Parity = 'none';
esp_port.Terminator = 'CR/LF';
smcAxis =1;

load('PowerCalibration800nm_20170712.mat');
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);

fclose(esp_port);
fopen(esp_port);



%deg = rot_calibr(0.5,f1);smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(deg));fprintf(laser,smcErrCmnd); % set absolute position

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.05,f1))));

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.1,f1))));


fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.2,f1))));

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.3,f1))));

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(0.5,f1))));

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(1.0,f1))));

fprintf(esp_port,strcat(num2str(smcAxis),'PA',num2str(rot_calibr(1.5,f1))));



fprintf(esp_port,'1PA0');


fclose(esp_port);
delete(esp_port);
