% This programe runs through the mat files generated previously and check
% if the heart rate generation is good

% first sample is great

% second sample was recorded continuesly
% middle part of 2.2 is not great but should still be usable
% 2.3 is not great either
% several frames in 2.4 is totally not usable, but all of a sudden, it
% becomes very good. Is there a problem with the code?

% last part of 3.3 has problem, the else are pretty good

% I'm guessing if continues imaging is causing trouble. but need
% confirmation when processed image. 

clear;
load('C:\Users\Zhengyi\OneDrive\Zebrafish heart beat profile\3.4 extractedProfileFromMaxIntensityRegion.mat')

for n = 1:size(profile,1)
    plot(profile(n,:));
    pause(0.1);
end