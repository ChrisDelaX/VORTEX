clear all
close all
da=0.01;
aExt=0:da:pi;
aInt=pi:-da:0;
zoom = 10;
myaxis = [-1 zoom -1 zoom];
newFig
for r=.5:2*zoom
    x=[(r+.25).*cos(aExt) (r-.25).*cos(aInt)];
    y=[(r+.25).*sin(aExt) (r-.25).*sin(aInt)];
    patch(x,y,'k')
    patch(x,-y,'k')
end
axis equal
axis(myaxis)
fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv2.eps', '-r300');
