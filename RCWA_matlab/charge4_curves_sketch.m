clear all
close all

Lb_max = 1.42;%1;%
alph = 3;%   %slope in degrees
depth = 5.6;% %grating depth
F_max = 0.51;
Lb_min = 2*depth*tan(deg2rad(alph))/(1-F_max);
%Lb_min=Lb_max*1e-4;%sqrt(.5);%0.8;

F_moy=.45;%.48;%.1;%
Lb_rat = Lb_min/Lb_max;%.5;%0.9;%.5;%   % ratio Lb_min/Lb_max

r_max=40;%10;%5;%
nb_pts=20;  % for one period Lb


% Figures
% -------
figure
hold on
%grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.35,0.35,0.35];
%mygrey = [.55,.55,.55];
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
set(gca,'visible','off')
axis equal

tic
depart_tmp=now;
depart=datestr(depart_tmp)

i=0;
j=0;
Lb_step = Lb_max;
r_ext=Lb_step*(1+F_moy)/4;
r_int=Lb_step*(1-F_moy)/4;
sin_max=1;
sin_min=Lb_rat^(1/2);
while 2*r_int*sin_min < r_max
    
    while 2*r_int*sin_min < r_max
        
        a_min = 2*asin(sin_min);

        %r_ext
        if 2*r_ext*sin_max <= r_max
            a_max = 2*asin(sin_max);
        elseif 2*r_ext*sin_min < r_max
            a_max = 2*asin(r_max/(2*r_ext));
        else
            a_max = a_min;
            r_ext=r_max/(2*sin_min);
        end
        nb_th=ceil((a_max-a_min)*r_ext/Lb_max*nb_pts+1);
        th_ext=linspace(a_min,a_max,nb_th);
        r=r_ext*ones(1,nb_th);
        [x,y] = pol2cart(th_ext,r);
        x=-x+r_ext;

        %r_int
        a_max = 2*asin(sin_max);
        if 2*r_int*sin_max <= r_max
            if 2*r_ext*sin_max > r_max
                x(end+1)=r_max*sin_max;
                y(end+1)=r_max*cos(asin(sin_max));
            end
            a_max = 2*asin(sin_max);
        else
            a_max = 2*asin(r_max/(2*r_int));
        end
        nb_th=ceil((a_max-a_min)*r_int/Lb_max*nb_pts+1);
        th_int=linspace(a_max,a_min,nb_th);
        r=r_int*ones(1,nb_th);
        [x2,y2] = pol2cart(th_int,r);
        x2=-x2+r_int;

        x(end+1:end+nb_th) = x2;
        y(end+1:end+nb_th) = y2;
        patch(x,y,mygrey,'edgecolor',mygrey);
        patch(-x,y,mygrey,'edgecolor',mygrey);
        patch(x,-y,mygrey,'edgecolor',mygrey);
        patch(-x,-y,mygrey,'edgecolor',mygrey);
%         if j==0
%             patch(x,-y,mygrey,'edgecolor',mygrey);
%         end
        i=i+1;
        r_ext=Lb_step*(1+F_moy+2*i)/4;
        r_int=Lb_step*(1-F_moy+2*i)/4;
%         if i == 3
%             break
%         end
    end    
    %increments
    i=0;
    j=j+1;
    Lb_step = Lb_max/Lb_rat^(j);
    r_ext=Lb_step*(1+F_moy+2*i)/4;
    r_int=Lb_step*(1-F_moy+2*i)/4;
    sin_max=Lb_rat^(j/2);
    sin_min=Lb_rat^((j+1)/2);
%     if j == 2
%         break
%     end
end
nb_dc=j-1  % number of discontinuities


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
str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',Lb_max,depth,Lb_rat,F_moy);
title(str,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')
%tick2latex

%save
print('-depsc2',sprintf('curves_sketch_Lbrat=%3.2f_F=%3.2f_r=%d.eps',Lb_rat,F_moy,r_max), '-r300');

% %Full screen
% pos0=get(0,'Screensize');
% pos1=get(gcf,'Position');
% pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
% set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]










