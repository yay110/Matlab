clear;
load('C:\Users\Zhengyi\OneDrive\Zebrafish heart beat profile\3.3b extractedProfileFromMaxIntensityRegion.mat')


%define parameters of the signal for FFT
Fs = 20;                        %sampling frequency
T = 1/Fs;                       % sampling period
L = sectionLength*Fs;       % Length of signal
t = (0:L-1)*T;                  %Time vector
profile = profile0;

for n =1:size(profile,1)
    subProfile = profile(n,:);
    Y = fft(subProfile);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    peaks = findpeaks(P1);
    beatFrequency(n) = f(find(P1==max(peaks)));
    normP1 = P1/max(peaks);
    
    %                  plot(f,normP1)
end

subplot(1,2,1)
plot(beatFrequency)
ylim([0,5])

subplot(1,2,2)
plot(profile(1,:))