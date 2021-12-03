x = 0:50:5000;

figure;
title('dev g = - 10%')
ylabel('signal changes in portion to B0')

ylabel('\DeltaS')
xlabel('imposed b-values (s/mm^2)')
hold on

y = exp(-x*0.95^2*0.0007)-exp(-x*0.0007);
plot(x,y,'LineWidth',2)
y = exp(-x*0.81*0.0007)-exp(-x*0.0007);
plot(x,y,'LineWidth',2)
y = exp(-x*0.85^2*0.0007)-exp(-x*0.0007);
plot(x,y,'LineWidth',2)

title([])
legend({'\DeltaL -5%','\DeltaL -10%','\DeltaL -15%'})
box on
set(gca,'LineWidth',2)
set(gca,'Fontsize',14)
set(gca,'yTick',[0 0.06 0.12])
