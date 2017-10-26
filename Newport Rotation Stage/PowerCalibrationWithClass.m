%% Setting up initial parameters for the com port
stage = NewportESP();
laser = Chameleon('COM1');
laser.Wavelength = 1050;
laser.SOpen;
%% Power Calibration
npoints = 30;
power = zeros(npoints,1);
rotation = power;
max_rotation = 90;
step = floor(max_rotation/npoints);
for i=1:npoints
    rotation(i) = (i-1)*step;
    stage.absPosition = rotation(i);
    power(i) = input('Please input the power reading: ');
end

stage.absPosition= 0;
%%make sure all the points are right before saving
today = datevec(date);
fname = sprintf('doubleWP-PowerCalibration1050nm_%d%02d%02d.mat',today(1),today(2),today(3));
save(fname,'rotation','power');
figure;
plot(rotation,power,'o');
ft = fittype('A*sin(2*pi()*x/T +phi) + offset','independent',{'x'}, 'coefficients',{'A','T','phi','offset'});
f1 = fit(rotation,power,ft,'StartPoint',[0.8 80 -10 1]);
hold on;
plot(f1)

%%

stage.Close;
%status = esp_port.Status;
laser.Close;
