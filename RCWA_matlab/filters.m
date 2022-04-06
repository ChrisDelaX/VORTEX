%clear all;close all;
%warning off MATLAB:singularMatrix
%warning off MATLAB:break_outside_of_loop

load filters.mat

datestr(now)


close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width

% 
band = [3250 4250]/1000;%[11 13.2];%
bandtick = [3250 3500 3750 4000 4250]/1000;%[11.0 11.5 12.0 12.5 13.0];%
% liss_lb_t = [lb_t(1):1e-4:lb_t(end)];


% Figures
% -------
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
fnz = 'Arial'; % fontname
fsz = 26; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2.2;  % line width 
%bandx = [lb_min lb_max];


% Filter WIDE
% -----------

dessin_filter
grid on
x1 = filter_wide(:,1);
y1 = filter_wide(:,2);
%f1 = plot(x1,y1,'k','linewidth',lwz);
p_wide = fit(x1,y1,'smoothingspline');
x2 = 3.25:.01:4.25;
y2 = p_wide(x2);
f2 = plot(x2,y2,'color',mygreen,'linewidth',2*lwz);
leg = legend('  filter wide');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
tit = title(['     Filter WIDE']);
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','HorizontalAlignment','center')

print('-depsc2',sprintf('filter_wide.eps'), '-r300');


% Filter WIDE+3800
% ----------------

dessin_filter
x1 = filter_3800(:,1);
y1 = filter_3800(:,2);%
p_3800 = fit(x1,y1,'smoothingspline');
p_wide_3800 = fit(x1,y1.*p_wide(x1),'smoothingspline');
x2 = 3.25:.01:4.25;
plot(x2,p_wide(x2),'--','color',myred,'linewidth',lwz);
plot(x2,p_3800(x2),'k','linewidth',lwz);
plot(x2,p_wide_3800(x2),'color',mygreen,'linewidth',2*lwz);
leg = legend('  Wide filter','  3800 filter','  Wide + 3800');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
set(leg,'box','on','linewidth',lwz)
tit = title(['      WIDE + 3800']);
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','HorizontalAlignment','center')

print('-depsc2',sprintf('filter_3800.eps'), '-r300');

%break


% Filter WIDE+4040
% ----------------

dessin_filter
x1 = filter_4040(:,1);
y1 = filter_4040(:,2);%.*p_wide(x1);
p_4040 = fit(x1,y1,'smoothingspline');
p_wide_4040 = fit(x1,y1.*p_wide(x1),'smoothingspline');
x2 = 3.25:.01:4.25;
plot(x2,p_wide(x2),'--','color',myred,'linewidth',lwz);
plot(x2,p_4040(x2),'k','linewidth',lwz);
plot(x2,p_wide_4040(x2),'color',mygreen,'linewidth',2*lwz);
leg = legend('  Wide filter','  4040 filter','  Wide + 4040');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
set(leg,'box','on','linewidth',lwz)
tit = title(['     WIDE + 4040']);
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','HorizontalAlignment','center')

print('-depsc2',sprintf('filter_4040.eps'), '-r300');

