clear all
close all
warning off MATLAB:polyfit

figure('name','ge')
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
h=1000e-6;
l=2:.01:30;
bandx = [0 30];
load data_hawkins.mat


% coeff n
%--------
subplot(2,2,1)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[3.8 4.2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Refractive index    n','FontSize',fsz)

T=300;
A=-6.040e-3*T + 11.05128 ; B= 9.295e-3*T + 4.00536 ; C = -5.392e-4*T + 0.599034 ; D = 4.151e-4*T + 0.09145; E = 1.51408*T + 3426.5 ;
nn=(A+(B.*(l.^2)./((l.^2)-C))+(D.*(l.^2)./((l.^2)-E)));
n1=sqrt(nn);
%n1(1030:end)=3.41565:(3.415-3.41565)/(size(n1,2)-1030):3.415;
plot(l,n1,'color',[1 .2 0],'linewidth',lwz)
T=150;
A=-6.040e-3*T + 11.05128 ; B= 9.295e-3*T + 4.00536 ; C = -5.392e-4*T + 0.599034 ; D = 4.151e-4*T + 0.09145; E = 1.51408*T + 3426.5 ;
nn=(A+(B.*(l.^2)./((l.^2)-C))+(D.*(l.^2)./((l.^2)-E)));
n2=sqrt(nn);
%n2(1040:end)=3.3917:(3.391-3.3917)/(size(n2,2)-1040):3.391;
plot(l,n2,'color',[0 .5 0],'linewidth',lwz)
T=50;
A=-6.040e-3*T + 11.05128 ; B= 9.295e-3*T + 4.00536 ; C = -5.392e-4*T + 0.599034 ; D = 4.151e-4*T + 0.09145; E = 1.51408*T + 3426.5 ;
nn=(A+(B.*(l.^2)./((l.^2)-C))+(D.*(l.^2)./((l.^2)-E)));
n3=sqrt(nn);
%n3(1040:end)=3.3757:(3.375-3.3757)/(size(n3,2)-1040):3.375;
plot(l,n3,'color',[0 .2 1],'linewidth',lwz)
leg=legend(' 300 K',' 150 K','  50 K');
set(leg,'Fontname',fnz,'FontSize',0.7*fsz)%,'location','southwest')
set(leg,'box','on','linewidth',lwz)

R1=((n1-1)./(n1+1)).^2;
T1=1-R1;
R2=((n2-1)./(n2+1)).^2;
T2=1-R2;
R3=((n3-1)./(n3+1)).^2;
T3=1-R3;

bandx = [10 30];
ll=12:.01:30;

% coeff k
%--------
subplot(2,2,2)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
%set(gca,'ylim',[3*1e-6 2*1e-3])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Extinction coefficient    K','FontSize',fsz)

k1=pchip(data_ge_k(:,1),data_ge_k(:,2),ll);
plot(ll,k1,'color',[1 .2 0],'linewidth',lwz)
k2=pchip(data_ge_k(:,1),data_ge_k(:,5),ll);
plot(ll,k2,'color',[0 .5 0],'linewidth',lwz)
k3=pchip(data_ge_k(:,1),data_ge_k(:,7),ll);
plot(ll,k3,'color',[0 .2 1],'linewidth',lwz)
leg=legend(' 300 K',' 150 K','  50 K');
set(leg,'Fontname',fnz,'FontSize',0.7*fsz,'location','northwest')
set(leg,'box','on','linewidth',lwz)

alpha=4*pi.*k1./(ll*1e-6);
abs1=ones(1,size(l,2));
abs1(size(l,2)-size(k1,2)+1:end)=exp(-alpha.*h);
alpha=4*pi.*k2./(ll*1e-6);
abs2=ones(1,size(l,2));
abs2(size(l,2)-size(k2,2)+1:end)=exp(-alpha.*h);
alpha=4*pi.*k3./(ll*1e-6);
abs3=ones(1,size(l,2));
abs3(size(l,2)-size(k3,2)+1:end)=exp(-alpha.*h);


% Reflectance
%-----------
subplot(2,2,3)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[30 60])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Reflectance  (%)','FontSize',fsz)

mRabs1=(1-R1).*abs1;
Rabs1=R1.*abs1;
Rtot1 = R1 .* (1 + (mRabs1.^2)./(1-Rabs1.^2));
plot(ll,Rtot1(size(l,2)-size(ll,2)+1:end).*100,'color',[1 .2 0],'linewidth',lwz)
mRabs2=(1-R2).*abs2;
Rabs2=R2.*abs2;
Rtot2 = R2 .* (1 + (mRabs2.^2)./(1-Rabs2.^2));
plot(ll,Rtot2(size(l,2)-size(ll,2)+1:end).*100,'color',[0 .5 0],'linewidth',lwz)
mRabs3=(1-R3).*abs3;
Rabs3=R3.*abs3;
Rtot3 = R3 .* (1 + (mRabs3.^2)./(1-Rabs3.^2));
plot(ll,Rtot3(size(l,2)-size(ll,2)+1:end).*100,'color',[0 .2 1],'linewidth',lwz)
leg=legend(' 300 K',' 150 K','  50 K');
set(leg,'Fontname',fnz,'FontSize',0.7*fsz,'location','southwest')
set(leg,'box','on','linewidth',lwz)




% Transmittance
%--------------
subplot(2,2,4)
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
%set(gca,'ylim',[20 60])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.7*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Transmittance  (%)','FontSize',fsz)

Ttot1 = (1-R1).*mRabs1 ./ (1-Rabs1.^2);
plot(ll,Ttot1(size(l,2)-size(ll,2)+1:end).*100,'color',[1 .2 0],'linewidth',lwz)
Ttot2 = (1-R2).*mRabs2 ./ (1-Rabs2.^2);
plot(ll,Ttot2(size(l,2)-size(ll,2)+1:end).*100,'color',[0 .5 0],'linewidth',lwz)
Ttot3 = (1-R3).*mRabs3 ./ (1-Rabs3.^2);
plot(ll,Ttot3(size(l,2)-size(ll,2)+1:end).*100,'color',[0 .2 1],'linewidth',lwz)
leg=legend(' 300 K',' 150 K','  50 K');
set(leg,'Fontname',fnz,'FontSize',0.7*fsz,'location','southwest')
set(leg,'box','on','linewidth',lwz)


print('-depsc2',sprintf('TH_ge.eps'), '-r300');


break 



