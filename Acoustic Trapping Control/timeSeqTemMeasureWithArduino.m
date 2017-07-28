%connect arduino
clear;

board = Arduino();
board.connect;

time = 300;
temperature(1:time,1:2) = 0;

for t = 1:time
    t
    tic
    temperature(t,:) = board.temperature;
% temperature(t,:)
timeSpend = toc
    pause(0-toc);
end

plot(temperature(:,1),'r');
hold on;
plot(temperature(:,2),'b');
% plot(temperature(:,3),'g');
% plot(temperature(:,4),'k');
% plot(temperature(:,5),'c');
hold off

ylim([20 70])

board.disconnect;