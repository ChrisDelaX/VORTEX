clear all
close all

Lb_max = 1.42;%1;%
alph = 3;%   %slope in degrees
depth = 5.6;% %grating depth
F_moy=.48;%.45;%.1;%

r_max=10;%30;%
nb_pts=100;  % for one period Lb


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
%set(gcf, 'color', 'none','inverthardcopy', 'off');
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(0,'DefaultTextInterpreter', 'latex')

tic
depart_tmp=now;
depart=datestr(depart_tmp)

i=0;
j=0;
r_ext=Lb_max*(.5+F_moy/2);
r_int=Lb_max*(.5-F_moy/2);





while r_ext < r_max
    nb_th=ceil(r_ext/Lb_max*nb_pts+1);
    th_ext=linspace(0,2*pi,nb_th);
    r=r_ext*ones(1,nb_th);
    [x,y] = pol2cart(th_ext,r);    
    nb_th=ceil(r_int/Lb_max*nb_pts+1);
    th_int=linspace(0,2*pi,nb_th);
    r=r_int*ones(1,nb_th);
    [x2,y2] = pol2cart(th_int,r);
    x(end+1:end+nb_th) = x2;
    y(end+1:end+nb_th) = y2;
    patch(x,y,mygrey,'edgecolor',mygrey);    
    i=i+1;
    r_ext=Lb_max*(.5+F_moy/2+i);
    r_int=Lb_max*(.5-F_moy/2+i);
end




%depart
%fin
toc
%datestr(now)


% Graphics
la_tmp=.05;%
line([-r_max r_max]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
text(r_max*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-r_max*(1+5*la_tmp),0,'$180^{\circ}$')
line([-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),[-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(r_max/sqrt(2)*(1+la_tmp),r_max/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-r_max/sqrt(2)*(1+4*la_tmp),-r_max/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
line([0 0],[-r_max r_max]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-r_max*la_tmp,r_max*(1+2*la_tmp),'$90^{\circ}$')
text(-r_max*la_tmp,-r_max*(1+2*la_tmp),'$270^{\circ}$')
line([r_max/sqrt(2) -r_max/sqrt(2)]*(1+la_tmp),[-r_max/sqrt(2) r_max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-r_max/sqrt(2)*(1+4*la_tmp),r_max/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(r_max/sqrt(2)*(1+la_tmp),-r_max/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
axis equal
xlabel('$\mu m$')
ylabel('$\mu m$')
set(gca,'xLim',[-r_max r_max]*(1+8*la_tmp))
set(gca,'ylim',[-r_max r_max]*(1+6*la_tmp))
tick2latex

%save
print('-depsc2',sprintf('curves_agpm_F=%3.2f_r=%d.eps',F_moy,r_max), '-r300')
%print('-dpng',sprintf('curves_agpm_F=%3.2f_r=%d.png',F_moy,r_max), '-r300')


%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]










