global T0 R0 absor

%clear all
%close all



% data
%-----
h=1000e-6;
%l=.18:.01:4.75;
bandxx = [0 6];
infrasil_data;
ll=[.18 ,.23  ,.29  ,.4  ,.54 ,.8    ,1   ,1.2 ,1.5 ,1.9 ,2.2 ,2.5  ,3   ,3.4 ,3.9   ,4.7]';
kk=[3e-5,1e-6,1e-7,7e-8,8e-8,1.4e-7,8e-8,4e-8,7e-8,4e-7,8e-7,1.5e-6,4e-6,8e-6,5e-5,4e-4]';
%plot(ll,kk)


% coeff n
%--------
n=pchip(infra_l,infra_n,l);
R0=((n-1)./(n+1)).^2;
T0=1-R0;

% coeff k
%--------
k=pchip(ll,kk,l);
alpha=4*pi.*k./(l*1e-6);
absor=exp(-alpha.*h);

% Reflectance
%-----------
mRabs=(1-R0).*absor;
Rabs=R0.*absor;
Rtot = R0 .* (1 + (mRabs.^2)./(1-Rabs.^2));

% Transmittance
%--------------
Ttot = (1-R0).*mRabs ./ (1-Rabs.^2);




% FIGURES
% -------

% figure('name','sio2')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   1000   800],'PaperUnits','points','PaperPosition',[400    0   1000   800]); % paper 
% fnz = 'Arial'; % fontname
% fsz = 20; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width

% subplot(2,2,1)
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandxx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[1 1.6])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
% ylabel('Refractive index    n','FontSize',fsz)
% plot(l,n,'color',[0 .5 0],'linewidth',lwz)

% subplot(2,2,2)
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandxx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[1e-8 10])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
% ylabel('Extinction coefficient    K','FontSize',fsz)
% plot(l,k,'color',[0 .5 0],'linewidth',lwz)

% subplot(2,2,3)
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandxx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[0 10])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
% ylabel('Reflectance  (%)','FontSize',fsz)
% plot(l,Rtot.*100,'color',[0 .5 0],'linewidth',lwz)

% subplot(2,2,4)
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandxx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[90 100])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
% ylabel('Transmittance  (%)','FontSize',fsz)
% plot(l,Ttot.*100,'color',[0 .5 0],'linewidth',lwz)


%print('-depsc2',sprintf('TH_sio2.eps'), '-r300');
