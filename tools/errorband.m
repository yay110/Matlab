% errorBand
% y = 0;
x = y(:,[1,3,4, 5,6,7]);
control = 0;

x(x==0) = nan;

meanX = nanmean(x, 2);
stdX = nanstd(x,0, 2);

t = (1:length(meanX))';

fill([t;flipud(t)],[meanX-stdX;flipud(meanX+stdX)],[.8 .8 .8],'linestyle','none');
hold on;
plot(t,meanX,'LineWidth',2,'Color','k');
ylim([0 1.5 ]);
xlim([0 length(t)]);

label = [20, 80];
for i = 1:length(label);    
plot([label(i) label(i)], [0 2], 'k');
end

% control = control/mean(control);
% plot(control);
for i =1:size(x,2)

s = plot(x(:,i),'k');
   s.Color(4) = 0.5;
end

hold off;
xlabel('Time (minutes)');
ylabel('Beat Rate Change (\Deltaf/f)');
title('Verapamil Treatment');
