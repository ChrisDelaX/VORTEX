clear all
close all
warning off MATLAB:polyfit

figure('name','diam')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   1000   800],'PaperUnits','points','PaperPosition',[400    0   1000   800]); % paper 
fnz = 'Arial'; % fontname
fsz = 20; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2.2;  % line width


% data
%-----
h=300e-6;
l=.55:.01:20;
%l=4.4:.1:5.2;
bandx = [0 20];
ll=[.55 , .75 ,1    ,1.5 ,2   ,2.25,2.45,2.56,2.68,2.72,2.81,2.93,3.01,3.08,3.67,3.70,3.85,3.95,4.065,4.10,4.23,4.44,4.61,4.7 ,4.72,4.76,4.9,4.94,4.98,5.05,5.2,5.25,5.4 ,5.5 ,5.7,6   ,7   ,8   ,13  ,18  ,20]';
%tt=[.710,.7105,.711,.713,.715,.705,.657,.657,.666,.668,.66 ,.657,.689,.687,.593,.5765,.578,.58 ,.628,.46 ,.4  ,.458,.464,.463,.4 ,.393,.394,.399,.48,.5  ,.565,.601,.65,.676,.689,.692,.701,.707,.709]';
tt=[.691,.6988,.7035,.708,.711,.713,.713,.705,.657,.657,.666,.668,.66 ,.657,.689,.687,.593,.5765,.578,.58 ,.628,.46 ,.4  ,.458,.464,.463,.4 ,.393,.394,.399,.48,.5  ,.565,.601,.65,.676,.689,.692,.701,.707,.709]';


% coeff n
%--------
subplot(2,2,1)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[2.37 2.425])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Refractive index    n','FontSize',fsz)

lb=(l.*1e3).^2;
A=1 ; B=0.3306 ; C=175.0 ; D= 4.3356 ; E=106.0 ; 
nn=(A + (B.*lb./(lb-C^2)) + (D.*lb./(lb-E^2)));
n=sqrt(nn);
plot(l,n,'color',[0 .5 0],'linewidth',lwz)
R=((n-1)./(n+1)).^2;
T=1-R;


% Transmittance
%--------------
subplot(2,2,4)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[35 80])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Transmittance  (%)','FontSize',fsz)

Ttot=pchip(ll,tt,l);
plot(l,Ttot.*100,'color',[0 .5 0],'linewidth',lwz)


% coeff k
%--------
subplot(2,2,2)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
%set(gca,'ylim',[90 100])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Extinction coefficient    K','FontSize',fsz)

lnn = -T.^2 + sqrt(T.^4 + 4.*R.^2.*Ttot.^2);
lnd = 2.*Ttot.*R.^2;
alpha=-log(lnn./lnd)./h;
k=l.*1e-6.*alpha./(4*h);
plot(l,k,'color',[0 .5 0],'linewidth',lwz)


% Reflectance
%-----------
subplot(2,2,3)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[0 45])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Reflectance  (%)','FontSize',fsz)

abs=exp(-alpha.*h);
mRabs=(1-R).*abs;
Rabs=R.*abs;
Rtot = R .* (1 + (mRabs.^2)./(1-Rabs.^2));
plot(l,Rtot*100,'color',[0 .5 0],'linewidth',lwz)




print('-depsc2',sprintf('TH_diam.eps'), '-r300');
