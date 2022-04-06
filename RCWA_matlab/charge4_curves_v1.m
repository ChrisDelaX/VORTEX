clear all
close all

nb_r=10;%30;%
%nb_theta = 14;%28;%51;%    %per quadrant
scale=.1;%.05;

F=0.45;%.1;%

lp=4;
nb_r_0 = nb_r;
r_max=nb_r;%30;%12;%

r = linspace(0,1,nb_r);
% thetmp = linspace(0,pi/2,nb_theta)';
% 
% theta = sin(thetmp)*pi/2;
% theta(end+1:2*nb_theta)=pi-theta(end:-1:1);
% theta(end+1:4*nb_theta)=2*pi-theta(end:-1:1)


% Figures
% -------
figure
hold on
%grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.4,0.4,0.4];
fnz = 'Arial'; % fontname: Helvetica
fsz = 12; % fontsize: 10
fwz = 'normal';  % fontweight: bold
msz = 8; % marker size
lwz = 2;  % line width
set(gcf,'color',mywhite)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(0,'DefaultTextInterpreter', 'latex')


theta=2*pi*(0:99)'/99;

tic
depart_tmp=now;
depart=datestr(depart_tmp)
for i=nb_r:-1:1
    r=(2*i-1+F)/4*ones(100,1);
    [x,y] = pol2cart(theta,r);
    x=x-(2*i-1+F)/4;
    %x=x-i/2;
    patch(x,y,mygrey,'edgecolor',mygrey);
    patch(-x,y,mygrey,'edgecolor',mygrey);
    r=(2*i-1-F)/4*ones(100,1);
    [x,y] = pol2cart(theta,r);
    x=x-(2*i-1-F)/4;
    %x=x-i/2+F/2;
    patch(x,y,'w','edgecolor','w');
    patch(-x,y,'w','edgecolor','w');
end

%depart
%fin
toc
%datestr(now)

% Graphics
la_tmp=.05;%
line([-r_max r_max]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
text(r_max*(1+1.5*la_tmp),0,'$0^{\circ}$')%,'Interpreter','latex')
text(-r_max*(1+5*la_tmp),0,'$180^{\circ}$')%,'Interpreter','latex')
line([-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),[-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(r_max/sqrt(2)*(1+la_tmp),r_max/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')%,'Interpreter','latex')
text(-r_max/sqrt(2)*(1+4*la_tmp),-r_max/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')%,'Interpreter','latex')
line([0 0],[-r_max r_max]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-r_max*la_tmp,r_max*(1+2*la_tmp),'$90^{\circ}$')%,'Interpreter','latex')
text(-r_max*la_tmp,-r_max*(1+2*la_tmp),'$270^{\circ}$')%,'Interpreter','latex')
line([r_max/sqrt(2) -r_max/sqrt(2)]*(1+la_tmp),[-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-r_max/sqrt(2)*(1+4*la_tmp),r_max/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')%,'Interpreter','latex')
text(r_max/sqrt(2)*(1+la_tmp),-r_max/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')%,'Interpreter','latex')
axis equal
xlabel('$\mu m$')%,'Interpreter','latex')
ylabel('$\mu m$')%,'Interpreter','latex')
%set(gca,'xLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'ylim',[-r_max r_max]*(1+6*la_tmp))%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
%str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',Lb_max,depth,Lb_rat,F_moy);
%tit=title(str);%,'Interpreter','latex');
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
plotTickLatex2D%('xlabeldy',100)

%save
print('-depsc2',sprintf('curves_v1_r=%d.eps',nb_r), '-r300');


%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]

