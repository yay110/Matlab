close all
clear
clc


%% Input
Clipoff             = 50;
Dist                = [300 350 400 500 600 700]; %Distance from fiber output in [mm]

%%
for k = 1:length(Dist)
    fn                  = ['Martin' num2str(Dist(k)/10) '_' num2str(Clipoff) '.xls'];
    [num,txt,raw]       = xlsread(fn);
    t                   = num(:,1);
    Pow(k,:)            = [mean(num(:,2)) std(num(:,2))]*1e-3;
    Pos.X(k,:)          = [mean(num(:,3)) std(num(:,3))]*1e-3;
    Pos.Y(k,:)          = [mean(num(:,4)) std(num(:,4))]*1e-3;
    Wid.V(k,:)          = [mean(num(:,5)) std(num(:,5))]*1e-3;
    Wid.W(k,:)          = [mean(num(:,6)) std(num(:,6))]*1e-3;
end
d                   = 0:1000;
cV                  = fit(Dist',Wid.V(:,1),'poly1');
Wid.cV              = cV.p1*d+cV.p2;
cW                  = fit(Dist',Wid.W(:,1),'poly1');
Wid.cW              = cW.p1*d+cW.p2;
w                   = mean([cV.p2 cW.p2])*sqrt(2)/sqrt(log(2));
slope               = mean([cV.p1 cW.p1])*sqrt(2)/sqrt(log(2));

%% lenses
L = [25.4, 50:25:150, 200:50:300, 400, 500:250:1000];
m = 0;
for l = 1:length(L)
    for k = 1:length(L)
        if L(k)>L(l)
            m = m+1;
            M(m) = L(k)/L(l);
            Lk(m) = L(k);
            Ll(m) = L(l);
        end
    end
end

%% Plotting
% BeamPosSTD          = std(sqrt(Pos.X'.^2+Pos.Y'.^2));
figure(1)
    subplot(121),plot(Dist, Wid.V(:,1)*sqrt(2)/sqrt(log(2)),'*'),hold on,plot(Wid.cV*sqrt(2)/sqrt(log(2)),'r'),xlabel('[mm]'),legend Data Fit
    subplot(122),plot(Dist, Wid.W(:,1)*sqrt(2)/sqrt(log(2)),'*'),hold on,plot(Wid.cW*sqrt(2)/sqrt(log(2)),'r'),xlabel('[mm]'),legend Data Fit
BeamWidMEAN.V             = Wid.V(:,1);
BeamWidMEAN.W             = Wid.W(:,1);

figure(2)
plot(Dist, Wid.V(:,1)*sqrt(2)/sqrt(log(2)),'*'),hold on,plot(Wid.cV*sqrt(2)/sqrt(log(2)),'r'),xlabel('[mm]'),legend Data Fit
plot(Dist, Wid.W(:,1)*sqrt(2)/sqrt(log(2)),'*'),hold on,plot(Wid.cW*sqrt(2)/sqrt(log(2)),'r'),xlabel('[mm]'),legend Data Fit
    plot(d,slope*d+w),xlabel('Distance from output [mm]'),ylabel('\omega_0 [mm]')



