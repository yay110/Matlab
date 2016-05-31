% Rapid acquisition


%% freqency and amplitude for the function generator
signalFrequency = 20;
signalAmplitude = 5;
signalOffset = 2.5;
phaseOffset = 0;

%% initialisation

% camera initialisation
% vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode');
%
% triggerconfig(vid,'manual')
% src = getselectedsource(vid);
% vid.FramesPerTrigger=1;
% vid.TriggerRepeat=Inf;
%
% start(vid)
% pause(1)
%
% src.ExposureTime = 0.01; %in seconds
%
% frameNumber = 50;
% pixels = 2048;
% pixelH = pixels;
% pixelV = pixels;
%
% frames = uint16(zeros(pixelH,pixelV,frameNumber));
% times=zeros(1,frameNumber);
%
% disp('Hamamatsu camera initialisation complete');

% function generator initialisation
fungen = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x05E6::0x3390::1242585::0::INSTR', 'Tag', '');
if isempty(fungen)
    fungen = visa('NI', 'USB0::0x05E6::0x3390::1242585::0::INSTR');
else
    fclose(fungen);
    fungen = fungen(1);
end

fopen(fungen);

fprintf(fungen, ['APPL:SIN ', num2str(signalFrequency), 'HZ, ', num2str(signalAmplitude), 'VPP, ',num2str(signalOffset),'V']);
fprintf(fungen, ['BURST:MODE TRIGGERED;NCYCLES INF;PHASE ',num2str(phaseOffset)]);
fprintf(fungen, 'TRIGGER:SOURCE EXT');

disp('Function Generator initialisation complete');

% Tunable lens initialisation
lens = Optotune('COM9');
lens.Open;
lens.setModeFrequency(signalFrequency);
lens.setModeLowerCurrent(0);
lens.setModeUpperCurrent(lens.max_current);





fprintf(fungen, 'BURST:STATE ON');
lens.sinusoidalMode;

% refresh the burst mode to eliminate the clock difference between function
% generator and the tunable lens. 
stop = 0;
while ~stop,
    delay(5);
    fprintf(fungen, 'BURST:STATE OFF');
    fprintf(fungen, ['APPL:SIN ', num2str(signalFrequency), 'HZ, ', num2str(signalAmplitude), 'VPP, ',num2str(signalOffset),'V']);
    fprintf(fungen, ['BURST:MODE TRIGGERED;NCYCLES INF;PHASE ',num2str(phaseOffset)]);
    fprintf(fungen, 'TRIGGER:SOURCE EXT');
    fprintf(fungen, 'BURST:STATE ON');
end


%% close the connections
% delete(vid);
% fclose(fungen);
% lens.Close;