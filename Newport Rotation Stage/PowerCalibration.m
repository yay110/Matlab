%% Setting up initial parameters for the com port
esp_port = serial('COM11');
esp_port.Baudrate=19200;
%esp_port.StopBits=1;
esp_port.Parity='none';
esp_port.Terminator='CR/LF';

%% Open port
fopen(esp_port);
status = esp_port.Status;

%% Choose initial settings
smcAxis =1;
%
smcErrCmnd = strcat(num2str(smcAxis),'ID?');
smcErr = query(esp_port,smcErrCmnd);

%%
smcErrCmnd = '2TP'; smcErr = query(esp_port,smcErrCmnd); %% check absolute position
% smcErrCmnd = '1VA8';smcErr = query(esp_port,smcErrCmnd); % set movement
% spee

smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(0.0)); fprintf(esp_port,smcErrCmnd); % set absolute position
%fprintf(obj,'cmd')

%fwrite(esp_port, 'tb');
%pause(2);
%bts = esp_port.BytesAvailable;
%response = fread(esp_port);

%% Power Calibration
npoints = 30;
power = zeros(npoints,1);
rotation = power;
max_rotation = 90;
step = floor(max_rotation/npoints);
for i=1:npoints
    rotation(i) = (i-1)*step;
    smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(rotation(i))); fprintf(esp_port,smcErrCmnd); % set absolute position
    power(i) = input('Please input the power reading: ');
end

fprintf(esp_port,'1PA0');
%%make sure all the points are right before saving
today = datevec(date);
fname = sprintf('PowerCalibration850nm_%d%02d%02d.mat',today(1),today(2),today(3));
save(fname,'rotation','power');
figure;
plot(rotation,power,'o');
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);
hold on;
plot(f1)

%%

fclose(esp_port);
%status = esp_port.Status;

