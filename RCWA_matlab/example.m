clear all
close all
figure
hold on
grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
set(gcf,'color',mywhite)
%set(gca,'box','on','linewidth',2)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(0,'DefaultTextInterpreter','latex')
xlabel('This is my Xlabel.')
ylabel('This is my Ylabel.')

x=9:-1:4;
y=-1:4;
plot(x,y,'color',mygreen,'linewidth',2)
%set(gca,'xaxislocation','top')   %optionnel
%set(gca,'yaxislocation','right')   %optionnel
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')

tick2latex

%save
print('-depsc2',sprintf('example.eps'), '-r300');

