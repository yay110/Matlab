%
%
%
function convolutionExampleFT()
    tRange=[-5:0.1:4.9]; dt=diff(tRange(1:2));
    f=@(t) (t==0 | t==2);
    g=@(t) exp(-t.^2./(2*0.5^2));
    
    % calculate the transfer function
    G=fft(ifftshift(g(tRange)));
    % convolve by multiplying
    F=fft(f(tRange));
    FxG=F.*G;
    fxg=ifft(FxG);
    
    % for display, calculate the frequency after Fourier transform
    df=1./(dt*numel(tRange)); % the highest will be determined by the Nyquist theorem
    fRange=df*([1:numel(tRange)]-1-floor(numel(tRange)./2));
    
    %
    % Output
    %
    close all;
    figure('Position',[100 100 1024 768]);
    subplot(3,2,1);
    plot(fRange,abs(fftshift(F))); title('F(f)');
    subplot(3,2,3);
    plot(fRange,abs(fftshift(G))); title('G(f)');
    subplot(3,2,5);
    plot(fRange,abs(fftshift(FxG))); title('F.G(f)');
    
    subplot(3,2,2);
    plot(tRange,f(tRange)); title('f(t)');
    subplot(3,2,4);
    plot(tRange,g(tRange)); title('g(t)');
    subplot(3,2,6);
    plot(tRange,fxg); title('fxg(tau)');
end