%
%
%
function convolutionExamples()
    tRange=[-5:0.1:4.9];
    
    f=@(t) (t==0 | t==2);
    g=@(t) exp(-t.^2./(2*0.5^2));
    
    % calculate the convolution by integrating for one or more values of tau
    function result=fxg(tauRange)
        dt=diff(tRange(1:2));
        result=zeros(size(tauRange));
        for tauIdx=1:numel(tauRange),
            tau=tauRange(tauIdx);
            result(tauIdx)=sum(f(tRange).*g(tRange-tau))*dt;
        end
    end
    
    %
    % Output
    %
    close all;
    figure('Position',[100 100 1024 768]);
    axs(1)=subplot(3,1,1);
    plot(tRange,f(tRange)); title('f(t)');
    axs(2)=subplot(3,1,2);
    plot(tRange,g(tRange)); title('g(t)');
    subplot(3,1,3);
    axs(3)=plot(tRange,fxg(tRange)); title('fxg(tau)');
    
    linkaxes(axs)
end