

close all

figure('FileName','Diamant');
set(gcf,'color',[1 1 1])
set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images

fs = 16; % fontsize
band = [11 13.2];



transmittance = subplot(2,2,1:2);
hold on
grid on
set(gca,'FontSize',fs,'FontWeight','bold')
set(gca,'XLim',band)
set(gca,'ylim',[0.4 1],'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorGrid','on')

tmp12 = zeros(size(lb_t)) + mean(tmp1+tmp2)/2;
TsansZOG=2./(n2lb+1);
plot(lb_t,(tmp1+tmp2)./2,'k',lb_t,tmp12,'k--',lb_t,TsansZOG,'k-.','Linewidth',2)
xlabel('Wavelength \lambda  (µm)')
ylabel('Transmittance')
l=legend(' T (\lambda)',' mean T',' no ZOG','location','northeast');
set(l,'FontSize',fs*.7)





retardance = subplot(2,2,3:4);
hold on
grid on
set(gca,'YMinorTick','on')
set(gca,'YMinorGrid','on')

[AX,H1,H2] = plotyy(lb_t,retard,lb_t,retard*360);
retmoy = zeros(size(lb_t)) + mean(retard);
plot(lb_t,retmoy,'k-.','Linewidth',2)
%plot(lb_t,zeros(size(lb_t)) + .5,'k--')
%axis tight
set(AX(1),'FontSize',fs,'FontWeight','bold','XLim',band,'XMinorTick','on','YColor','k')
set(AX(2),'FontSize',fs,'FontWeight','bold','XLim',band,'XMinorTick','on','YColor','k','ylim',get(AX(1),'ylim')*360)
set(AX(1),'ylim',[0.47 0.53],'ytick',[.47 .48 .49 .5 .51 .52 .53])
set(AX(2),'ylim',get(AX(1),'ylim').*360,'ytick',[170 175 180 185 190])
set(H1,'color','k','linewidth',2)%,'marker','.')
set(H2,'color','k','linewidth',2)%,'marker','.')
xlabel('Wavelength \lambda  (µm)')
ylabel(AX(1),'Retardance (Waves)')
ylabel(AX(2),'Phase Shift (°)')
l=legend(' \Delta\Phi_{TE-TM} (\lambda)',' mean \Delta\Phi_{TE-TM}','location','northeast');
set(l,'FontSize',fs*.7)


% pserror = subplot(2,4,5:6);
% hold on
% grid on
% set(gca,'FontSize',fs,'FontWeight','bold')
% set(gca,'XLim',band)
% set(gca,'XMinorTick','on')
% set(gca,'YScale','log')
% 
% plot(lb_t,smooth(abs(tmp3-pi)),'k.-','Linewidth',2)
% errmoy = zeros(size(lb_t)) + mean(abs(tmp3-pi));
% plot(lb_t,errmoy,'k--','Linewidth',1)
% xlabel('Wavelength \lambda  (µm)')
% ylabel('absolute error  | \Delta\Phi_{TE-TM}|  (rad)')



saveas(gcf,['spectrum_RMS_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.fig']);
print('-dpng', ['spectrum_RMS_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.png'], '-r300');

