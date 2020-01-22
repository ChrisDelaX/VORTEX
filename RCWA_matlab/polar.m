clear all
close all
newFig
set(gca,'visible','off')
axis equal
lp=2;
maskstr='SGVC';
r=5;
x=0;
y=0;
ang=0:0.001:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'k--','linewidth',lwz)
quiver(0,0,1,0,'color',myblue,'linewidth',lwz,'MaxHeadSize',1.2)
quiver(0,0,-1,0,'color',myblue,'linewidth',lwz,'MaxHeadSize',0)
nTH=16;
TH=0:2*pi/nTH:2*pi;
for i=1:nTH
    x0=r*cos(TH(i));
    y0=r*sin(TH(i));
    x1=(cos(lp*TH(i)));
    y1=(sin(lp*TH(i)));
    quiver(x0,y0,x1,y1,'color',myred,'linewidth',lwz,'MaxHeadSize',1.2)
    quiver(x0,y0,-x1,-y1,'color',myred,'linewidth',lwz,'MaxHeadSize',0)
    x2=(cos(lp/2*TH(i)));
    y2=(sin(lp/2*TH(i)));    
    quiver(x0,y0,x2,y2,'k--','linewidth',lwz/2,'MaxHeadSize',0)
    quiver(x0,y0,-x2,-y2,'k--','linewidth',lwz/2,'MaxHeadSize',0)
end
text(0,2*r/9,'input','FontSize',1.5*fsz,'HorizontalAlignment','center')
text(0,r/9,'polarization','FontSize',1.5*fsz,'HorizontalAlignment','center')

text(0,12*r/9,'{\bf Charge-2 vortex}','FontSize',1.5*fsz,'HorizontalAlignment','center')
%text(0,12*r/9,'{\bf Charge-4 vortex}','FontSize',1.5*fsz,'HorizontalAlignment','center')
%text(0,12*r/9,'{\bf Charge-6 vortex}','FontSize',1.5*fsz,'HorizontalAlignment','center')


print('-depsc2',sprintf('polar_%s%d.eps',maskstr,lp), '-r300');
