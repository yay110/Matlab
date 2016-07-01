% % Rapid acquisition
% function [] = fastACQ(nrCycles)
%
% if nargin <1
%     nrCycles = 1;
% end

%% freqency and amplitude for the function generator
generatorVoltage = [-0.127,0.168];  % Voltage on the function generator, unit Volt.
phaseOffset = 96;                   % Offset between tunable lens and the scanning mirror.
% Vary with frequency

signalFrequency = 10;
currentETL = [50,150];

signalAmplitude = generatorVoltage(2)-generatorVoltage(1);
signalOffset = mean(generatorVoltage);

pixels = 256;
%% initialisation

% camera initialisation
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');

triggerconfig(vid,'hardware','risingEdge','EdgeTrigger');
src = getselectedsource(vid);


vid.ROIPosition = [(2048-pixels)/2 (2048-pixels)/2 pixels pixels];
src.ExposureTime = 0.01/2048*pixels; %in seconds

% each trigger will include 2 full cycle, i.e. four volumes of images
vid.FramesPerTrigger=1/signalFrequency/src.ExposureTime*2;
%vid.FramesPerTrigger = 10;
%vid.TriggerRepeat=Inf;

% nrFrames = vid.FramePerTrigger;
% pixelH = pixels;
% pixelV = pixels;
%
% frames = uint16(zeros(pixelH,pixelV,nrFrames));

disp('Camera Ready');

%% function generator initialisation
fungen = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x05E6::0x3390::1242585::0::INSTR', 'Tag', '');
if isempty(fungen)
    fungen = visa('NI', 'USB0::0x05E6::0x3390::1242585::0::INSTR');
else
    fclose(fungen);
    fungen = fungen(1);
end

fopen(fungen);

if signalFrequency == 1
    fprintf(fungen, ['APPL:RAMP ', num2str(signalFrequency), 'HZ, ', num2str(signalAmplitude), 'VPP, ',num2str(signalOffset),'V']);
    fprintf(fungen, 'FUNCTION:RAMP:SYMM 50');
else
    fprintf(fungen, ['APPL:SIN ', num2str(signalFrequency), 'HZ, ', num2str(signalAmplitude), 'VPP, ',num2str(signalOffset),'V']);
end
fprintf(fungen, ['BURST:MODE TRIGGERED;NCYCLES 100;PHASE ',num2str(phaseOffset)]);
fprintf(fungen, 'TRIGGER:SOURCE EXT');

disp('Function Generator ready');

%% Tunable lens initialisation
lens = Optotune('COM9');
lens.Open;

if signalFrequency == 1
    lens.triangularMode;
else
    lens.sinusoidalMode;
end
lens.setModeFrequency(signalFrequency);
lens.setModeLowerCurrent(currentETL(1));
lens.setModeUpperCurrent(currentETL(2));
disp('Tunable lens ready');

fprintf(fungen, 'BURST:STATE ON');


pause(15);
start(vid)
disp('Start Acquiring Frames!');
pause(5)

frames = getdata(vid);

%% close the connections
delete(vid);
%fprintf(fungen, 'BURST:STATE OFF');
fclose(fungen);
lens.Close;
disp('All devices disconnected succesfully')


%binning the frames by factor of binning and display


outputFileName = 'C:\acquisitionData\img_stack2.tif';
binning = 1;
newPixels = pixels/binning;
newFrames = zeros(newPixels,newPixels,vid.FramesPerTrigger,'uint16');
for n = 1:vid.FramesPerTrigger
    c = frames(:,:,:,n);
    %      c = reshape(c,[binning newPixels binning newPixels]);
    %      c = sum(c,1);
    %      c = sum(c,3);
    %      c = reshape(c,[newPixels newPixels]);
    %     c = c/max(max(c))*65535;
    newFrames(:,:,n) = c;
    imwrite(c, outputFileName, 'WriteMode', 'append');
end

% video display
%implay(newFrames);

