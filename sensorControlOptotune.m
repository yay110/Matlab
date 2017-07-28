% lens = Optotune('COM9');
% lens = lens.Open();
%
% analogInput = lens.getAnalogInput();
%
% % Lookup table;
% index = 200 / lens.max_bin;

nLoop = 100;

times = zeros(1,nLoop);
for n=1:nLoop
    tic;
%     lens.getAnalogInput();
    current = lens.analogInput*200/lens.max_bin;
%     lens.setCurrent(current);
    times(n) = toc;
end

figure;

plot(times*1e3);
% lens = lens.Close();
