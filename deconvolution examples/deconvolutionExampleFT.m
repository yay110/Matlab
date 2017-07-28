%
%
%
function deconvolutionExampleFT()
    tRange=[-5:0.1:4.9]; dt=diff(tRange(1:2));
    f=@(t) (t==0 | t==2);
    g=@(t) exp(-t.^2./(2*0.5^2));
    
    % calculate the transfer function
    G=fft(ifftshift(g(tRange)));
    % convolve by multiplying
    F=fft(f(tRange));
    FxG=F.*G;
    fxg=ifft(FxG);
    clear F FxG; % in reality we don't know these
    
    % determine a filter
    NSR=10/100;
    Hw=conj(G)./(abs(G).^2+NSR.^2);
    
    % deconvolve
    FxG=fft(fxg);
    F_dec=FxG.*Hw;
    f_dec=ifft(F_dec);
    
    % for display, calculate the frequency after Fourier transform
    df=1./(dt*numel(tRange)); % the highest will be determined by the Nyquist theorem
    fRange=df*([1:numel(tRange)]-1-floor(numel(tRange)./2));
    
    %
    % Output
    %
    close all;
    figure('Position',[100 100 1024 768]);
    subplot(3,2,1);
    plot(fRange,abs(fftshift(FxG))); title('FxG(f)');
    subplot(3,2,3);
    plot(fRange,abs(fftshift(Hw))); title('Hw(f)');
    subplot(3,2,5);
    plot(fRange,abs(fftshift(F_dec))); title('F_{dec}(f)');
    
    subplot(3,2,2);
    plot(tRange,fxg); title('fxg(t)');
    subplot(3,2,4);
    plot(tRange,real(fftshift(ifft(Hw)))); title('hw(t)');
    subplot(3,2,6);
    plot(tRange,f_dec); title('f_{dec}(t)');
end