clear all
close all
da=0.01;
width = .25;%.1;%
aExt = pi/4:da:pi;
aInt = pi:-da:pi/4;
zoom = 10;
myaxis = [-1 zoom -1 zoom];

newFig
for r=.25:.5/1.22:10
    x=[(r+width).*cos(aExt)-r*sqrt(2) (r-width).*cos(aInt)-r*sqrt(2)];
    y=[(r+width).*sin(aExt) (r-width).*sin(aInt)];
    patch([-1-width*sqrt(2) -1+width*sqrt(2) zoom+width*sqrt(2) zoom-width*sqrt(2)],[-1 -1 zoom zoom],'k')
    patch([8 zoom zoom 8],[8 8 zoom zoom],'k')
    patch(x,y,'k');patch(y,x,'k')
    patch(x,-y,'k');patch(-y,x,'k')
    patch(-x,y,'k');patch(y,-x,'k')
    patch(-x,-y,'k');patch(-y,-x,'k')
end
axis equal
axis(myaxis)
fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv6a.eps', '-r300');

%break

newFig
i=0;
for r=.25:.5:17
%    rr=(1+i/2)*r;
%    rr=r*tan(pi/2/12*i);
    rr=r/cos(pi/2/16*i);
    x=[(rr+width).*cos(aExt)-rr*sqrt(2) (rr-width).*cos(aInt)-rr*sqrt(2)];
    y=[(rr+width).*sin(aExt) (rr-width).*sin(aInt)];
    patch([-1-width*sqrt(2) -1+width*sqrt(2) zoom+width*sqrt(2) zoom-width*sqrt(2)],[-1 -1 zoom zoom],'k')
    patch([1-width*sqrt(2) 1+width*sqrt(2) -1 -1],[-1 -1 1+width*sqrt(2) 1-width*sqrt(2)],'k')
    patch(x,y,'k');patch(y,x,'k')
    patch(x,-y,'k');patch(-y,x,'k')
    patch(-x,y,'k');patch(y,-x,'k')
    patch(-x,-y,'k');patch(-y,-x,'k')
    i=i+1;
end
axis equal
axis(myaxis)
fontsize=16;
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv6b.eps', '-r300');

break
newFig
for r=.25:.5:50
    x=[(r+width/2).*cos(aExt)-(r+width/2)*sqrt(2) (r-width/2).*cos(aInt)-(r-width/2)*sqrt(2)];
    y=[(r+width/2).*sin(aExt) (r-width/2).*sin(aInt)];
    %patch([-1 0 1],[zoom 1 zoom],'k')
    patch(x,y,'k');patch(y,x,'k')
    patch(x,-y,'k');patch(-y,x,'k')
    patch(-x,y,'k');patch(y,-x,'k')
    patch(-x,-y,'k');patch(-y,-x,'k')
end
axis equal
axis(myaxis)
set(gca,'fontsize',fontsize,'tickdir','out')
tick2latex;print('-depsc2','sgv6c.eps', '-r300');



