% Multi-param (F,d)
%------------------

close all

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
lwz = 2;  % line width


figure('FileName','Diamant');
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper
hold on
%grid on
set(gca,'box','on','linewidth',2)

% grid
% x=0.35:0.05:0.55;
% y=3.5:0.5:7;
% y=y';
% z=zeros(size(x,2),size(y,1));
% hG = surfc(x,y,z');

%hS = pcolor(Xparam,Yparam,log10(fmulti'));
hS = surfc(Xparam,Yparam,log10(fmulti'));
shading interp
set(hS,'EdgeColor', 'none')
%shading faceted
%set(hS,'edgecolor','k')
hC = colorbar;
set(hC,'box','on','linewidth',2)
%set(hC,'location','northoutside')
%set(hC,'YScale','log')
%set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
%colormap bone% gray%
contourf(peaks(100))

%xlabel('Average Filling Factor F','Fontname',fnz,'FontWeight',fwz,'FontSize',fsz*1.2)
xlabel('Filling Factor F','FontSize',12)
ylabel('Grating Depth h (\mum)','FontSize',12)
set(gca,'YDir','reverse')
%set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')

xmin = Xparam(1);
xmax = Xparam(end);
ymin = Yparam(1);
ymax = Yparam(end);
axis([xmin xmax ymin ymax])

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

if zoom(1) == 'o'
    set(gca,'xlim',[.348 .572])
    set(gca,'ylim',[3.47 7.03])
elseif zoom(1) == 'i'
    set(gca,'xlim',[.348 .512])
    set(gca,'ylim',[3.47 6.08])
end
%set(gca,'xtick',[.35 .4 .45 .5 .55])
%set(gca,'ytick',[4 4.25 4.5 4.75 5 5.25])
