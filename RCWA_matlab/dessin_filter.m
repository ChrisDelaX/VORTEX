%close all

figure('name','filter')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',band,'xtick',bandtick)
%set(gca,'YScale','log')
set(gca,'ylim',[0 1.2])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('Transmission','FontSize',fsz)


