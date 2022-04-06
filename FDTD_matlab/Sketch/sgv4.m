clear all
close all
da=0.01;
width = .25;%.1;%
aExt = 0:da:pi;
aInt = pi:-da:0;
zoom = 10;
myaxis = [-1 zoom -1 zoom];

newFig
for r=.25:.5:20
    x=[(r+width).*cos(aExt)-r (r-width).*cos(aInt)-r];
    y=[(r+width).*sin(aExt) (r-width).*sin(aInt)];
    patch([-width -width width width],[-1 zoom zoom -1],'k')
    patch([-1 -1 zoom 1],[2 zoom zoom 2],'k')
    patch(x,y,'k')
    patch(x,-y,'k')
    patch(-x,y,'k')
    patch(-x,-y,'k')
end
axis equal
axis(myaxis)
fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv4a.eps', '-r300');


newFig
i=0;
for r=.25:1:8
%    rr=(1+i/2)*r;
%    rr=r*tan(pi/2/10*i);
    rr=r/cos(pi/2/8*i);
    x=[(rr+width).*cos(aExt)-rr (rr-width).*cos(aInt)-rr];
    y=[(rr+width).*sin(aExt) (rr-width).*sin(aInt)];
    patch([-width -width width width],[-1 zoom zoom -1],'k')
    patch(x,y,'k')
    patch(x,-y,'k')
    patch(-x,y,'k')
    patch(-x,-y,'k')
    i=i+1;
end
axis equal
axis(myaxis)
fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv4b.eps', '-r300');

break
newFig
for r=.25:.5:50
    x=[(r+width/2).*cos(aExt)-(r+width/2) (r-width/2).*cos(aInt)-(r-width/2)];
    y=[(r+width/2).*sin(aExt) (r-width/2).*sin(aInt)];
    patch([-1 0 1],[zoom 1 zoom],'k')
    patch(x,y,'k')
    patch(x,-y,'k')
    patch(-x,y,'k')
    patch(-x,-y,'k')
end
axis equal
axis(myaxis)
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv4c.eps', '-r300');



