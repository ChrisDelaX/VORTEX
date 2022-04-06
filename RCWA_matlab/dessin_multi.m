% Multi-param (F,d)
%------------------

close all

figure
hold on
%grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.4,0.4,0.4];
fnz = 'Arial'; % fontname: Helvetica
fsz = 18; % fontsize: 10
fwz = 'normal';  % fontweight: bold
msz = 8; % marker size
lwz = 2;  % line width
set(gcf,'color',mywhite)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')



hS = surfc(xparam.vals,yparam.vals,log10(nuldpt'));
%shading('interp')
shading('flat') %flat, faceted
set(hS,'EdgeColor', 'none')
hC = colorbar;
set(hC,'box','on','linewidth',lwz)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)

CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
colormap(CMRmap)
caxis([-3.5 0])

xlabel([xparam.name ' (' (xparam.units) ')'],'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz)
ylabel([yparam.name ' (' (yparam.units) ')'],'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz)
tit = title(sprintf('%s band [%3.2f-%3.2f %s], %s = %3.2f %s', band, ...
    lam.range(1), lam.range(2), lam.units, zparam.name, zparam.val, ...
    zparam.units));%,'horizontalalignment','left')
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)

%set(gca,'YDir','reverse')
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')

set(gca,'ticklength',-.9*get(gca,'ticklength'))

axis([xparam.range(1) xparam.range(2) yparam.range(1) yparam.range(2)]);

%%

% % lower bound
% xmin_low = xmin;
% ymin_low = xmin*Lb/tan(pente);
% xmax_low = xmax;
% ymax_low = xmax*Lb/tan(pente);
% [b]=line([xmin_low xmax_low],[ymin_low ymax_low]);
% set(b,'Linewidth',lwz,'color',[1 1 1]) 
% % upper bound
% xmin_up = xmin;
% ymin_up = (1-xmin)*Lb/tan(pente);
% xmax_up = xmax;
% ymax_up = (1-xmax)*Lb/tan(pente);
% [b]=line([xmin_up xmax_up],[ymin_up ymax_up]);
% set(b,'Linewidth',lwz,'color',[0 0 0])

% if zoom(1) == 'o'
%     set(gca,'xlim',[.348 .572])
%     set(gca,'ylim',[3.47 7.03])
% elseif zoom(1) == 'i'
%     set(gca,'xlim',[.348 .512])
%     set(gca,'ylim',[3.47 6.08])
% end

%set(gca,'xtick',[.35 .4 .45 .5 .55])
%set(gca,'ytick',[2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8])
%set(gca,'xtick',[ .4 .45 .5 ])
%set(gca,'ytick',[ 2.6 2.8 3.0 3.2 3.4 3.6])
%set(gca,'ytick',[ 2.6 2.8 3.0 3.2 3.4 3.6])

plain = -3.8:.2:-2.6;
dash = -2.4:.2:-2.0;
dashdots = -1.8:.2:-0;
[C01,C02] = contour(xparam.vals,yparam.vals,log10(nuldpt'), [plain], 'w-');
[C03,C04] = contour(xparam.vals,yparam.vals,log10(nuldpt'), [dash], 'w--');
[C05,C06] = contour(xparam.vals,yparam.vals,log10(nuldpt'), [dashdots], 'w-.');
clabel(C01,C02, 'Color','w','FontSize',fsz);
clabel(C03,C04, 'Color','w','FontSize',fsz);
clabel(C05,C06, 'Color','w','FontSize',fsz);


