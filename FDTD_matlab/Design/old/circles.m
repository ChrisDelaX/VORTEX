clear all
close all
period = 1.42;
F = .87/period;
size = 30;%10;
myaxis = [-size size -size size];

fig=figure;
set(gcf,'color','none')
set(gca,'box','off')
axis equal; axis(myaxis)

da=0.01;
aExt=0:da:pi;
aInt=pi:-da:0;
for r=.5:2*size/period
    x=[(r+F/2).*cos(aExt)*period (r-F/2).*cos(aInt)*period];
    y=[(r+F/2).*sin(aExt)*period (r-F/2).*sin(aInt)*period];
    patch(x,y,'r')
    patch(x,-y,'r')
end

fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
print('-depsc2','sgv2.eps', '-r1000');

%addpath('altmany-export_fig-d966721/')
%export_fig(fig,'sgv2','-png')
