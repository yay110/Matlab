close all
% clear
% Zhengyi Yang
% Use tracked results from BoundaryTracking.m to do PSD analsys and plot
% Required data : two arrays "c" and "r" contains the position information
%   for the bead
clc

%%
%  load('stiffness29Jan.mat')
power = [31.23, 27.40, 23.21, 18.06, 16.71, 21.72, 24.91, 35.30, 43.24, 62.08,37.69];

%% Constants
fps                     = 100;
p2m                     = 0.0833*1e-6;                      %Pixel to meter [m/pix]
T                       = 297.15;                           %Temperature [K] 24 degrees
kB                      = 1.3806488e-23;                    %Boltzmann constant [m^2*kg/(s^2*K)]
Radius                  = 2.5e-6;                        %Radius of bead = 2.5um;
Viscosity               = 0.892e-6;                      % Kinematic viscosity m2/s
Density                 = 997;                           % Density Kg/m3
StokesDrag              = 6*pi*Density*Viscosity*Radius;
Cons4Fitting            = kB*T/(pi^2*StokesDrag);

%%processing parameters
LowF         = 0.01;              %lowest freqency for fitting
HighF        = 10;               %highest frequency for fitting
PPoints2Proc = [1:4,6:10];      %5th set of data is excluded because doesn't look good
PPoints2Fit  = [1:4,6:9];       

for powerID=7%PPoints2Proc;%size(c,2)/10
    
    %%============== reshape the matrix ==================
    zPart                 = c(:,powerID*10-9:powerID*10);
    Z(powerID,:)          = reshape(zPart,1,size(zPart,1)*size(zPart,2));
    xPart                 = r(:,powerID*10-9:powerID*10);
    X(powerID,:)          = reshape(xPart,1,size(xPart,1)*size(xPart,2));
    
    %%============== stiffness by equipartition=================
    xStiffness(powerID)      	= 1e6*kB*T./var(X(powerID,:));                     %[pN/um]
    zStiffness(powerID)         = 1e6*kB*T./var(Z(powerID,:));                     %[pN/um]
    xStd(powerID)               = std(X(powerID,:));
    zStd(powerID)               = std(Z(powerID,:));
    %%============= FFT ==================================
    sx                  = X(powerID,:);
    [~,~,f,~,PSDx]      = transformFourier(fps,sx,2);
    
    sz                  = Z(powerID,:);
    [~,~,~,~,PSDz]      = transformFourier(fps,sz,2);
    
    %%=======least square fitting into Lorentzian-shaped power spectrum======
    [fitx,~,~]          = fitPSD(f,PSDx,Cons4Fitting,LowF,HighF);
    f0x                 = sqrt(coeffvalues(fitx));
    kappaPSDx(powerID)  = 2*pi*StokesDrag*f0x;
    ConfInt             = sqrt(confint(fitx));
    ePSDx(powerID)      = f0x-ConfInt(1);
    
    [fitz,LowN,HighN]   = fitPSD(f,PSDz,Cons4Fitting,LowF,HighF);
    f0z                 = sqrt(coeffvalues(fitz));
    kappaPSDz(powerID)  = 2*pi*StokesDrag*f0z;
    ConfInt             = sqrt(confint(fitz));
    ePSDz(powerID)      = f0z-ConfInt(1);
    
    %%if plot?
    if  1
        %%=============== plot ===============================
        %size and postion of the plot
        scrnsz            	= 0.5.*get(0,'ScreenSize');
        nhist               = 50;
        h=figure('Position',scrnsz);
        
        ha = tight_subplot(2,3,[.09 .1],[.1 .05],[0.05 .01]);
        %         h=figure;
        
        
        if powerID == 7
            %limits for PSD plot
            XLIM                = [f(1),f(end)];
            YLIM                = [0.1*min([PSDz;PSDx]),10*max([PSDz;PSDx])];
            %Limits for time plot
            t                	= 0:1/fps:(length(sz)-1)/fps;          %Time vector [s]
            m2um                = 1e6;                                  %meter to micrometer
            wMax                = m2um*max([max(abs(X(powerID,:))),max(abs(Z(powerID,:)))]);
            wLim                = 0.8*[-wMax, wMax];
        end
        
        %=================Time plot=================
        %         subplot(2,3,1),
        axes(ha(1));
        plot(t(1:1:0.1*length(t)),m2um*X(powerID,1:1:0.1*length(t)),'b')
        %         title('X direction Positional Trace'),
        xlabel('Time [s]'),ylabel('X Position [\mum]'),xlim([t(1),0.1*t(end)]),ylim(wLim)
        set(gca,'Xtick',[0,10,20,30,40,50,60]);
        axis square;
        %         subplot(2,3,4),
        axes(ha(4));
        plot(t(1:1:0.1*length(t)),m2um*Z(powerID,1:1:0.1*length(t)),'b')
        %         title('Z direction Positional Trace'),
        xlabel('Time [s]'),ylabel('Z Position [\mum]'),xlim([t(1),0.1*t(end)]),ylim(wLim)
        set(gca,'Xtick',[0,10,20,30,40,50,60]);
        axis square;
        
        
        
        %==================PSD plot================
        %         subplot(2,3,3),
        axes(ha(3));
        loglog(f,1e12*PSDx,'b')
        xlim(XLIM),ylim(1e12*YLIM)
        %         title('Power spectrum density plot on X direction')
        hold;
        loglog(f(LowN:HighN),1e12*1/2*Cons4Fitting./(f0x.^2+f(LowN:HighN).^2),'r','LineWidth',3);
        xlabel('Frequency [Hz]');ylabel('PSD [\mum^2/Hz]')
        set(gca,'Xtick',[1e-2,1e-1,1e0,1e1]);
        
        axis square;
        
        hold off;
        
        %         subplot(2,3,6),
        axes(ha(6));
        loglog(f,1e12*PSDz)
        xlim(XLIM),ylim(1e12*YLIM)
        %         title('Power spectrum density plot on Z direction')
        hold;
        loglog(f(LowN:HighN),1e12*1/2*Cons4Fitting./(f0z.^2+f(LowN:HighN).^2),'r','LineWidth',3);
        xlabel('Frequency [Hz]');ylabel('PSD [\mum^2/Hz]')
        set(gca,'Xtick',[1e-2,1e-1,1e0,1e1]);
        axis square;
        
        hold off;
        
        
        %======================histgram plot==========================
        %         subplot(2,3,2),
        axes(ha(2));
        his=histfit(m2um*X(powerID,:),nhist);
        set(his(1),'facecolor','b'); set(his(2),'color','r','LineWidth',3)
        %         title(sprintf('Histogram, std = %1.2f um',m2um*xStd(powerID))),
        xlabel('X Position [\mum]'),ylabel('Counts [#]'),xlim(wLim)
        set(gca,'Xtick',-3:1:3,'Ytick',0:1000:5000);
        axis square;
        
        
        %         subplot(2,3,5),
        axes(ha(5));
        his=histfit(m2um*Z(powerID,:),nhist);
        set(his(1),'facecolor','b'); set(his(2),'color','r','LineWidth',3)
        %         title(sprintf('Histogram, std = %1.2f um',m2um*zStd(powerID))),
        xlabel('Z Position [\mum]'),ylabel('Counts [#]'),xlim(wLim)
        set(gca,'Xtick',-3:1:3,'Ytick',0:1000:5000);
        axis square;
        
        
        %  set(findall(gca,'type','text'),'FontSize',15,'fontWeight','bold')
        %save file
        %         saveas(h,num2str(powerID),'fig')
%         saveas(h,num2str(powerID),'pdf')
        %         %         plot2svg('myplot.svg',h,'png');
        %         saveas(h,num2str(powerID),'png')
        %         close all;
        
    end
end

AverKx=kappaPSDx(PPoints2Fit)./power(PPoints2Fit);
mean(AverKx)*1e6
std(AverKx)/sqrt(length(AverKx))*1e6
AverKz=kappaPSDz(PPoints2Fit)./power(PPoints2Fit);
mean(AverKz)*1e6
std(AverKz)/sqrt(length(AverKz))*1e6
%plot stiffness from equipartition
% figure
% subplot(2,2,1),plot(power, xStiffness,'o');xlim([0,max(power)*1.2]);ylim([0,max(xStiffness)*1.2])
% title('Stiffness on X direction from equipartition'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% subplot(2,2,2),plot(power, zStiffness,'o');xlim([0,max(power)*1.2]);ylim([0,max(zStiffness)*1.2])
% title('Stiffness on Z direction from equipartition'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% subplot(2,2,3),plot(1:10, xStiffness./power,'o');ylim([0,max(xStiffness./power)*1.1]);
% title('Stiffness on power on X direction from equipartition'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% subplot(2,2,4),plot(1:10, zStiffness./power,'o');ylim([0,max(zStiffness./power)*1.1]);
% title('Stiffness on power on Z direction from equipartition'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')

%plot stiffness from power spectrum density
% figure
% subplot(1,2,1),errorbar(power(PPoints2Proc), 1e6*kappaPSDx(PPoints2Proc),ePSDx(PPoints2Proc),'o');xlim([0,max(power(PPoints2Proc))*1.2]);ylim([0,max(1e6*kappaPSDx(PPoints2Proc))*1.2])
% title('Stiffness on X direction from PSD'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% hold on;
% %force to go through (0,0)
% slopex       = (power(PPoints2Fit)')\(kappaPSDx(PPoints2Fit)');
% plot(power(PPoints2Fit),1e6*power(PPoints2Fit)*slopex,'y','LineWidth',3);
% %poly fit
% fit=polyfit(power(PPoints2Fit),kappaPSDx(PPoints2Fit),1);
% plot(power(PPoints2Fit),1e6*(power(PPoints2Fit)*fit(1)+fit(2)),'r','LineWidth',3);
% hold off;
% 
% subplot(1,2,2),errorbar(power(PPoints2Proc), 1e6*kappaPSDz(PPoints2Proc),ePSDz(PPoints2Proc),'o');xlim([0,max(power(PPoints2Proc))*1.2]);ylim([0,max(1e6*kappaPSDz(PPoints2Proc))*1.2])
% title('Stiffness on Z direction from PSD'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% hold on;
% slopez       = (power(PPoints2Fit)')\(kappaPSDz(PPoints2Fit)');
% plot(power(PPoints2Fit),1e6*power(PPoints2Fit)*slopez,'r','LineWidth',3);
% hold off;
% 
% subplot(2,2,3),errorbar(1:10, 1e6*kappaPSDx./power,ePSDx./power,'o');ylim([0,max(1e6*kappaPSDx./power)*1.1]);
% title('Stiffness on power on X direction from PSD'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')
% subplot(2,2,4),errorbar(1:10, 1e6*kappaPSDz./power,ePSDz./power,'o');ylim([0,max(1e6*kappaPSDz./power)*1.1]);
% title('Stiffness on power on Z direction from PSD'),xlabel('Power [um]'),ylabel('Stiffness [pN/um]')