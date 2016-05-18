function[fit1,LowN,HighN]=fitPSD(f,PSD,Cons,LowF,HighF)
%Zhengyi Yang Feb. 2015
%least square fitting into Lorentzian-shaped power spectrum
%f          :frequency
%PSD        :Power spectrum denstiy
%LowF       :lowest freqency for fitting
%HighF      :Highest frequency for fitting

if nargin < 5
    HighF   = max(f);
end

if nargin < 4
    LowF    = min(f);
end

LowN        = find(f>LowF,1);
HighN       = find(f>HighF,1);

xdata   = f(LowN:HighN);
ydata   = PSD(LowN:HighN)/Cons;
ftype   = fittype('1/(a+x.^2)');
[fit1,~,fitinfo] = fit(xdata,ydata,ftype,'StartPoint',0);
%         residuals = fitinfo.residuals;
% a       = coeffvalues(fit1);
% f0      = sqrt(a);



