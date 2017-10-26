%% Setting up initial parameters for the com port
NewportESP = serial('COM5');
NewportESP.Baudrate=19200;
%esp_port.StopBits=1;
NewportESP.Parity='none';
NewportESP.Terminator='CR/LF';

%% Open port
fopen(NewportESP);
status = NewportESP.Status;

%% Choose initial settings
smcAxis =1;
%
smcErrCmnd = strcat(num2str(smcAxis),'ID?');
smcErr = query(NewportESP,smcErrCmnd);

%%
smcErrCmnd = '2TP'; smcErr = query(NewportESP,smcErrCmnd); %% check absolute position
% smcErrCmnd = '1VA8';smcErr = query(esp_port,smcErrCmnd); % set movement
% spee

smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(0.0)); fprintf(NewportESP,smcErrCmnd); % set absolute position
%fprintf(obj,'cmd')

%fwrite(esp_port, 'tb');
%pause(2);
%bts = esp_port.BytesAvailable;
%response = fread(esp_port);

%% Power Calibration
npoints = 30;
power = zeros(npoints,1);
rotation = power;
max_rotation = 60;
step = floor(max_rotation/npoints);
for i=1:npoints
    rotation(i) = (i-1)*step;
    smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(rotation(i))); fprintf(NewportESP,smcErrCmnd); % set absolute position
    power(i) = input('Please input the power reading: ');
end

fprintf(NewportESP,'1PA0');
%%make sure all the points are right before saving
today = datevec(date);
fname = sprintf('PowerCalibration750nm_%d%02d%02d.mat',today(1),today(2),today(3));
save(fname,'rotation','power');
figure;
plot(rotation,power,'o');
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);
hold on;
plot(f1)

%%

fclose(NewportESP);
%status = esp_port.Status;

