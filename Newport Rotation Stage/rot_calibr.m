%% Calibration
function deg = rot_calibr(targetPower,f1)
    deg = (asin((targetPower-f1.offset)/f1.A) - f1.phi)*f1.T/(2*pi());
    deg =deg-f1.T;
end