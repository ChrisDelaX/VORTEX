clear all
close all

nb_r = 15;%21;%
nb_theta = 16;%28;%51;%    %per quadrant
scale=.075;%.1;%

lp=4;
nb_r_0 = nb_r;

nb_theta = nb_theta+1;
r = linspace(0,1,nb_r);
thetmp = linspace(0,pi/2,nb_theta)';

theta = sin(thetmp)*pi/2;
theta(end+1:2*nb_theta)=pi-theta(end:-1:1);
theta(end+1:4*nb_theta)=2*pi-theta(end:-1:1);




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



tic
depart_tmp=now;
depart=datestr(depart_tmp)
for j=1:4*nb_theta
    j
    if abs(cos(theta(j)))>.1%.02%
        r = linspace(1/(2*nb_r+1),1,nb_r);
        for i=1:nb_r
            quiver(r(i)*cos(theta(j)),r(i)*sin(theta(j)),-r(i)*sin(lp/2*theta(j)),r(i)*cos(lp/2*theta(j)),scale/abs(cos(theta(j)))^(1/2),'k.')
            quiver(r(i)*cos(theta(j)),r(i)*sin(theta(j)),r(i)*sin(lp/2*theta(j)),-r(i)*cos(lp/2*theta(j)),scale/abs(cos(theta(j)))^(1/2),'k.')
            if i==1 && j==1
                fin=datestr(depart_tmp+toc/3600/24*nb_r*nb_theta)
            end
        end
    end
end

%depart
%fin
toc
%datestr(now)

% Graphics
axis equal
xlabel('$\mu m$')%,'Interpreter','latex')
ylabel('$\mu m$')%,'Interpreter','latex')
%set(gca,'xLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[-r_max r_max]*(1+6*la_tmp))%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
%str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',Lb_max,depth,Lb_rat,F_moy);
%tit=title(str);%,'Interpreter','latex');
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
plotTickLatex2D%('xlabeldy',100)

%save
print('-depsc2',sprintf('lines_v0_th=%d_r=%d.eps',nb_theta-1,nb_r), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



