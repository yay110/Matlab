            %% Setting up initial parameters for the com port
esp_port = serial('COM12');
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
smcErr = query(esp_port,smcErrCmnd)

%%
%smcErrCmnd = '1TP'; smcErr = query(esp_port,smcErrCmnd); %% check absolute position
% smcErrCmnd = '1VA8';smcErr = query(esp_port,smcErrCmnd); % set movement
% speed
%smcErrCmnd = '1PA1.0';smcErr = query(esp_port,smcErrCmnd); %% set absolute position
var = depth(m_index)/1000 + depth_low;
smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(var)); fprintf(esp_port,smcErrCmnd);

position_set = rot_calibr(p_max,f1); smcErrCmnd = strcat(num2str(smcAxis),'PA',num2str(position_set)); fprintf(esp_port,smcErrCmnd); % set absolute position
%fprintf(obj,'cmd')

%fwrite(esp_port, 'tb');
%pause(2);
%bts = esp_port.BytesAvailable;
%response = fread(esp_port);
%%

fclose(esp_port);
status = esp_port.Status;