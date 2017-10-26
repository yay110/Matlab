adata = cdata;

b = size(adata);
adata = adata';
adata(adata==0) = [];
adata = reshape(adata,b(2)-17*2,b(1));
adata = adata';

frequency = 20;
Fs = 2048/size(adata,1)*100; %sampling frequency;
T = 1/Fs;             % Sampling period       
L = size(adata,2);             % Length of signal
t = (0:L-1)*T;        % Time vector
f = Fs*(0:(L/2))/L;


for i = 1:size(adata,1)
    
X = adata(i,:);
Y = fft(X);

% P2 = abs(Y/L);
P2 = angle(Y);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

fftMap(i,:) = P1;
end

subplot(1,2,1)
imagesc(f,1:size(fftMap,1),fftMap)
% xticks(0:5:size(fftMap,2))

subplot(1,2,2)
plot(fftMap(:,find(f==20)));
