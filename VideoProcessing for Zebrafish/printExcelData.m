t1 = 0;
t2 = 0;
t3 = 0;
t4 = 0;


plot(t1(:,1).*24.*60,t1(:,2),'k-o' ,'MarkerFaceColor',[0 0 0],'MarkerSize',3)
ylim([1 3])
hold on
    

plot(t1(:,3).*24.*60,t1(:,4),'b-s', 'MarkerFaceColor',[0 0 1],'MarkerSize',5)



plot(t1(:,5).*24.*60,t1(:,6),'r-d', 'MarkerFaceColor',[1 0 0],'MarkerSize',8)

hold off;




