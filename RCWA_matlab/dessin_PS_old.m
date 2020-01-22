%lb_vect=[1.15 1.5 2 3.45 4.35 8.5 11];
lb_vect=[3.45 4.35 8.5 11];
npts_lb=size(lb_vect,2);
npts_E=61;%11;%3;%
pente_deg = pente/pi*180;

for band_nb=1:npts_lb    
    lb=lb_vect(band_nb);
    
%     for j=1:npts_E
%         E_man=E_man_min+(E_man_max-E_man_min)*(j-1)/(npts_E-1);
%         Lbzog = lb/E_man;


close all

%band_nb = 7;
[lb_vect(band_nb) lb_vect(band_nb)*1.2]

figure('FileName','Diamant');
set(gcf,'color',[1 1 1])
set(gcf,'Position', get(0,'Screensize')); % Maximize figure.
set(gcf,'PaperUnits','points','PaperPosition',[0 0 1440 878]) % screen sized images

nulldepth = subplot(2,3,1);
hold on
grid on
[b]=plot(Xparam(1,:),smooth(rms_err(band_nb,:),40),'Marker','.');
set(b,'Linewidth',2,'color',[1 0 0]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('rms error','FontSize',16,'FontWeight','bold')
%title('Null Depth','FontSize',16,'FontWeight','bold')
set(gca,'YScale','log')
set(gca,'FontSize',16,'FontWeight','bold')

phaseshift = subplot(2,3,4);
hold on
grid on
[b]=plot(Xparam(1,:),smooth(esp_phi(band_nb,:),40),'Marker','.');
set(b,'Linewidth',2,'color',[1 0 0]) 
% [c]=plot(Xparam(1,:),smooth(1./(nulldpt(2,:)*5e-3),4),'Marker','.');
% set(c,'Linewidth',2,'color',[0 1 0]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('Mean phase shift','FontSize',16,'FontWeight','bold')
%title('Null Depth','FontSize',16,'FontWeight','bold')
%set(gca,'YScale','log')
set(gca,'FontSize',16,'FontWeight','bold')

ParamY2 = subplot(2,3,2);
hold on
grid on
Yzog = lb_vect(band_nb)./Xparam(1,:);
plot(Xparam,smooth(Yzog,10),'--k')
[b]=plot(Xparam(1,:),Yparam2(band_nb,:),'Marker','.');
set(b,'Linewidth',2,'color',[0 0 1]) 
% Yzog = lb_vect(2)./Xparam(1,:);
% plot(Xparam,Yzog,'--k')
% [c]=plot(Xparam(1,:),Yparam2(2,:),'Marker','.');
% set(c,'Linewidth',2,'color',[0 1 0]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('Period (\mum)','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')

ParamY1 = subplot(2,3,3);
hold on
grid on
[b]=plot(Xparam(1,:),smooth(Yparam1(band_nb,:),2),'Marker','.');
set(b,'Linewidth',2,'color',[0 0 1]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('Filling factor','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')

ParamY3 = subplot(2,3,5);
hold on
grid on
[b]=plot(Xparam(1,:),smooth(Yparam3(band_nb,:),10),'Marker','.');
set(b,'Linewidth',2,'color',[0 0 1]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('Thickness (\mum)','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')

aspectRatio = subplot(2,3,6);
hold on
grid on
F_asp_rat = .5 - abs(.5-Yparam1(band_nb,:));
asprat = Yparam3(band_nb,:)./(F_asp_rat.*Yparam2(band_nb,:));
[b]=plot(Xparam(1,:),smooth(asprat,10),'Marker','.');
set(b,'Linewidth',2,'color',[0 0 1]) 
xlabel('Refraction Index','FontSize',16,'FontWeight','bold')
ylabel('Aspect Ratio','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')


saveas(gcf,['spectrum_pente_',sprintf('%3.2f',pente_deg),'_band_',sprintf('%3.2f',lb_vect(band_nb)),'-',sprintf('%3.2f',lb_vect(band_nb)*1.2),'.fig']);
print('-dpng', ['spectrum_pente_',sprintf('%3.2f',pente_deg),'_band_',sprintf('%3.2f',lb_vect(band_nb)),'-',sprintf('%3.2f',lb_vect(band_nb)*1.2),'.png'], '-r300');



    
end


