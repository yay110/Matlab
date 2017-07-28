%% This program will use the history power calibration data to predict power at new wavelength

wavelength = [720,800,850,900,960];
power = [0.85,1.76,1.51,1.16,0.6];

for i = 1:5
load(strcat('PowerCalibration',num2str(wavelength(i)),'nm_20170127.mat'));
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);

A(i)=f1.A;
T(i)=  f1.T;
phi(i) = f1.phi;
offset(i) = f1.offset;
end

power = [0.85,1.76,1.51,1.16,0.6];

% fitA = polyfit(power,A,1);
% plot(power,A,'r')
% hold;
% plot(power,fitA(1)*power+fitA(2),'b');
% hold;
% 
% fitoffset = polyfit(power,offset,1);
% plot(power,offset,'r')
% hold;
% plot(power,fitoffset(1)*power+fitoffset(2),'b');
% hold;

newPower = 1.76;
f2 = f1;
f2.T = mean(T);
f2.phi = mean(phi);
f2.A = fitA(1)*newPower+fitA(2);
f2.offset = fitoffset(1)*newPower+fitoffset(2);
f1 = f2;
