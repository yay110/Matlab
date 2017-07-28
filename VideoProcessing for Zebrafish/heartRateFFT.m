% clear;
mysnuggle=getenv('USERPROFILE');
% load(strcat(mysnuggle,'\OneDrive\Zebrafish heart beat profile\2.4c extractedProfileFromMaxIntensityRegion.mat'));
 profile = profile2;
 clear beatFrequency;


%define parameters of the signal for FFT
sectionLength = 9;
Fs = 20;                        %sampling frequency
T = 1/Fs;                       % sampling period
L = sectionLength*Fs;       % Length of signal
t = (0:L-1)*T;                  %Time vector

% reshape the matrix for 
for n=1:size(profile,1);
profile(n,:) = profile(n,:)/max(squeeze(profile(n,:)));
end
profile0=profile';
profile0=reshape(profile0,1,[]);
plot(profile0(1:180));
profile0(length(profile0)+1:ceil(length(profile0)/L)*L)=0;
% profile0= profile0(1:floor(length(profile0)/L)*L); 
profile0=reshape(profile0,L,[]);
profile = profile0';
clear profile0;


for n =1:size(profile,1)
    subProfile = profile(n,:);
    Y = fft(subProfile);
    P2 = abs(Y/L);
    P1(n,:) = P2(1:L/2+1);
    P1(n,2:end-1) = 2*P1(n,2:end-1);
    f = Fs*(0:(L/2))/L;
    peaks = findpeaks(P1(n,10:40));
    beatFrequency(n) = f(find(P1(n,:)==max(peaks)));
    normP1 = P1/max(peaks);
    
    %                  plot(f,normP1)
end
beatFrequency = beatFrequency';

subplot(1,2,1)
timeStamp = (1:size(beatFrequency,1));%*sectionLength/60;
plot(timeStamp,beatFrequency,'o')
% hold on
% beatFrequency0 = reshape(beatFrequency, [],3);
% beatFrequency0 = mean(beatFrequency0,2);
% timeStamp0 = (1:size(beatFrequency0,1));
% plot(timeStamp0,beatFrequency0,'*')
% hold off
ylim([0,6])

for n =75%1:size(beatFrequency,1)
%     n=116
subplot(2,2,2)
plot(profile(n,:))
subplot(2,2,4)
plot(f,P1(n,:))
pause(0.5)

end