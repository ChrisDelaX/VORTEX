

close all

figure('FileName','Diamant');
set(gcf,'color',[1 1 1])
set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images

fs = 16; % fontsize
band = [11 13.2];

nulldepth = subplot(2,4,1:2);
hold on
grid on
set(gca,'FontSize',fs,'FontWeight','bold')
set(gca,'XLim',band)
set(gca,'XMinorTick','on')
set(nulldepth,'YScale','log')
t=title(['                                                                                                                                           Perfomances of a diamond AGPM    [11-13.2 µm]    with    \Lambda = ',sprintf('%3.1f',Lb),' µm,    F = ',sprintf('%3.0f',F*100),' %,    h = ',sprintf('%3.1f',d),' µm    and   \alpha = ',sprintf('%3.2f',pente_deg),'°']);
set(t,'FontSize',20,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom')
%axis square
%gtext({['Slope \alpha: ',num2str(pente/pi*180),' °'],['Period \Lambda: ',num2str(Lb),' \mum'],['Depth h: ',num2str(d),' \mum'],['Filling factor F: ',num2str(F*100),' %']})

plot(lb_t,nul_res_sp_b,'k.-','Linewidth',2)
nulmoy = zeros(size(lb_t)) + null_res_sp;
plot(lb_t,nulmoy,'k--','Linewidth',1)
xlabel('Wavelength \lambda  (µm)')
ylabel('Null Depth')



transmittance = subplot(2,4,3:4);
hold on
grid on
set(gca,'FontSize',fs,'FontWeight','bold')
set(gca,'XLim',band)
set(gca,'XMinorTick','on')
set(gca,'YMinorGrid','on')

tmp12 = zeros(size(lb_t)) + mean(tmp1+tmp2)/2;
TsansZOG=2./(n2lb+1);
plot(lb_t,(tmp1+tmp2)./2,'k',lb_t,tmp12,'k--',lb_t,TsansZOG,'k-.','Linewidth',2)
xlabel('Wavelength \lambda  (µm)')
ylabel('Transmittance')
l=legend(' T(\lambda)',' mean T',' no ZOG','location','best');
set(l,'FontSize',fs*.7)


% reflectance = subplot(2,4,4);
% hold on
% grid on
% set(gca,'FontSize',fs,'FontWeight','bold')
% set(gca,'XLim',band)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorGrid','on')
% 
% tmp45 = zeros(size(lb_t)) + mean(tmp4+tmp5)/2;
% RsansZOG=(n2lb-1)./(n2lb+1);
% plot(lb_t,tmp4,'-',lb_t,tmp5,'--',lb_t,RsansZOG,'k-.','Linewidth',2)
% plot(lb_t,tmp45,'k--','Linewidth',1)
% xlabel('Wavelength \lambda  (µm)')
% ylabel('Reflectance')
% l=legend(' R_{TE}(\lambda)',' R_{TM}(\lambda)',' mean R',' no ZOG','location','best');
% set(l,'FontSize',fs*.7)


% dessinprofil = subplot(2,4,5);
% imagesc(prof)
% xlabel('Period \Lambda (%)')
% ylabel('Grating thickness (%)')


pserror = subplot(2,4,5:6);
hold on
grid on
set(gca,'FontSize',fs,'FontWeight','bold')
set(gca,'XLim',band)
set(gca,'XMinorTick','on')
set(gca,'YScale','log')

plot(lb_t,smooth(abs(tmp3-pi)),'k.-','Linewidth',2)
errmoy = zeros(size(lb_t)) + mean(abs(tmp3-pi));
plot(lb_t,errmoy,'k--','Linewidth',1)
xlabel('Wavelength \lambda  (µm)')
ylabel('absolute error  | \Delta\Phi_{TE-TM}|  (rad)')




retardance = subplot(2,4,7:8);
hold on
grid on
set(gca,'YMinorTick','on')
set(gca,'YMinorGrid','on')

[AX,H1,H2] = plotyy(lb_t,retard,lb_t,retard*360);
retmoy = zeros(size(lb_t)) + mean(retard);
plot(lb_t,retmoy,'k--','Linewidth',1)
%plot(lb_t,zeros(size(lb_t)) + .5,'k--')
%axis tight
set(AX(1),'FontSize',fs,'FontWeight','bold','XLim',band,'XMinorTick','on','YColor','k')
set(AX(2),'FontSize',fs,'FontWeight','bold','XLim',band,'XMinorTick','on','YColor','k','ylim',get(AX(1),'ylim')*360)
set(H2,'color','k','marker','.','linewidth',2)
xlabel('Wavelength \lambda  (µm)')
ylabel(AX(1),'Retardance (Waves)')
ylabel(AX(2),'Phase Shift (°)')



saveas(gcf,['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.fig']);
print('-dpng', ['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.png'], '-r300');

