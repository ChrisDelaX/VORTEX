%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Gestion de la mémoire
%---------------------
%---------------------

%Libération
%----------

save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
warning off MATLAB:polyfit


global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof retardghost dn_sp_T n_s_T n_p_T LFt

%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)
N=8;%12;%

%Allocation mémoire
%------------------

kx_mat=zeros(2*N+1);
ky_mat=zeros(2*N+1);
kIz_mat=zeros(2*N+1);
kIIIz_mat=zeros(2*N+1);
E=zeros(2*N+1);
A=zeros(2*N+1);
B=zeros(2*N+1);
D=zeros(2*N+1);
Omega=zeros(2*(2*N+1));
Tuu_l=zeros(2*(2*N+1));
Rud_l=zeros(2*(2*N+1));
Rdu_l=zeros(2*(2*N+1));
Tdd_l=zeros(2*(2*N+1));
F_lp1=zeros(4*(2*N+1));
X_lp1=zeros(2*(2*N+1));


%Lectures des entrées
%--------------------
%--------------------

%Choix du profil
%---------------

%1: profil rectangulaire simple
%   -> 1 couche (L=1)
%2: profil rectangulaire avec couche antireflet dans les creux et sur les bosses
%   -> 3 couches (L=3)
%3: profil rectangulaire avec couche antireflet sur les bosses uniquement
%   -> 2 couches (L=2)
%4: profil rectangulaire avec couche antireflet continue.5
%   -> 2 couches (L=2)
%5: profil rectangulaire simple avec couche d'arret continue
%   -> 2 couches (L=2)
%6: profil rectangulaire avec couche d'arret continue et antireflet dans les creux et sur les bosses
%   -> 4 couches (L=4)
%7: profil rectangulaire avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> 3 couches (L=3)
%8: profil rectangulaire avec couche d'arret continue et antireflet continue
%   -> 3 couches (L=3)
%9: profil rectangulaire simple avec couche d'arret discontinue
%   -> 2 couches (L=2)
%10: profil rectangulaire avec couche d'arret discontinue et antireflet dans les creux et sur les bosses
%   -> 4 couches (L=4)
%11: profil rectangulaire avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> 3 couches (L=3)
%12: profil rectangulaire avec couche d'arret discontinue et antireflet continue
%   -> 3 couches (L=3)

%13: profil trapézoidal simple
%   -> L couches
%14: profil trapézoidal avec couche antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2 couches
%15: profil trapézoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trapézoidal avec couche antireflet continue
%   -> L+1 couches
%17: profil trapézoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trapézoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trapézoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)

p=13;%1;%3;%

% %Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
% L=25;%50;%8;%
% pente_deg=0;%3.2;%2.75;% 
% pente=rad(pente_deg);

%pente_min = rad(4);
%pente_max = rad(7);
%pente=atan(0.1);%10/180*pi;   % => 10% = 5.71°
%pente=rad(0.57); % =1%
%pente=rad(0.86); % =1,5%
%pente=rad(1.15); % =2%;


%Choix matériau --> Permittivités
%--------------------------------

%1: Vide/Air
%2: CdTe
%3: Diamant
%4: Germanium
%5: Silicium
%6: ZnSe
%7: YF3
%8: Manuel
%9: AsGa
%10: n-laf32
%11: GASIR 2
%12: GASIR 1
%13: InP
%14: Infrasil
%15: KRS-5
%16: Si3N4
%17: ZnS
%18: n-lasf44
%19: ZnSe (T)
%20: LAH83

E_man=1.5;

%Permittivités milieu extérieur incident
EI_choice=1;

%Permittivités milieu extérieur émergent
EIII_choice=3;%5;%

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;%5;%

%Permittivités de la couche antireflet
E_AR_choice=1;
E_AR1_choice=1;
E_AR2_choice=1;

%Permittivités de la couche d'arret
E_arret_choice=1;


%Paramètres du réseau
%--------------------

%Pas du réseau (sublambda!)
% --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
Lb=.265;%
Lb_min=4.7;Lb_max=4.75;

%Epaisseurs µm
d=4.2;%4.3;%
d_min=12;%12.5;%
d_max=16;%14.5;%
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage
F=.36;%0.4;%
F_min=0.35;%0.42;%
F_max=0.45;%0.50;%

%Onde incidente
%--------------

%Longueur d'onde µm
lb=2.2;%
lb0=lb;
lb_min=11;%2;%.53;%.65;%1.15;%1.5;%
lb_max=13.2;%2.4;%.65;%.79;%1.4;%1.8;%
nlb=11;%21;%81;%

%Angle incidence non conique
theta=rad(0);
theta_min=rad(30); theta_max=rad(50);

%Angle incidence conique
phi=rad(0);
phi_min=rad(0); phi_max=rad(0);

%Angle polarisation
psi=rad(45);
psi_min=rad(0); psi_max=rad(0);




% Calcul rcwa
% -----------

F=0.4;%
Lb=4.6;
pente=rad(2.6);
L=16;
d=13.5489;%13.7315;%13.9565;%
nlb=81;
% nlb=11;%
% x0=[d];
% [x,fval,exitflag] = fminsearchbnd(@null,x0,d*.5,d*2,optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
% x
% [nuldpt,sigerr]=null(x)
% break

[nuldpt,sigerr]=null()



% Figures
% -------
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
fnz = 'Arial'; % fontname
fsz = 30; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2.2;  % line width 
bandx=[lb_min lb_max];
close all
fig0=figure('name','null')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   1400 600],'PaperUnits','points','PaperPosition',[400    0   2000 600]); % paper 
% trans
% ----------
fig3 = subplot(1,4,1);
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[0.8 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('TE-TM transmittance','FontSize',fsz)
plot(lb_t,tmp1,'--','color',myblue,'linewidth',lwz)
plot(lb_t,tmp2,'--','color',myred,'linewidth',lwz)
plot(lb_t,(tmp1+tmp2)./2,'k-','linewidth',lwz)
%plot(lb_t,(n2lb-1)./(n2lb+1),'k.:','linewidth',lwz)
leg=legend(' T_{TE}(\lambda)',' T_{TM}(\lambda)',' mean T');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','south')
set(leg,'box','on','linewidth',lwz)
tit=title('N-band AGPM  [11-13.2 µm]  with \Lambda = 4.6 µm,  F = 40%,  h = 13.5 µm  and  \alpha = 2.6°');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','horizontalalignment','left')
% refl
% -----------
fig4 = subplot(1,4,2);
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[0 0.2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('TE-TM reflectance','FontSize',fsz)
plot(lb_t,tmp4,'--','color',myblue,'linewidth',lwz)
plot(lb_t,tmp5,'--','color',myred,'linewidth',lwz)
plot(lb_t,(tmp4+tmp5)./2,'k-','linewidth',lwz)
leg=legend(' R_{TE}(\lambda)',' R_{TM}(\lambda)',' mean R');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
set(leg,'box','on','linewidth',lwz)
% Amp mismatch
% -----------
fig1 = subplot(1,4,3:4);
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[.9 1.05])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('TE-TM amplitude mismatch','FontSize',fsz)
plot(lb_t,tmp1./tmp2,'k-o','linewidth',lwz)
plot(lb_t,lb_t.*0+1,'k--','linewidth',lwz)
leg=legend(' intensity ratio  q(\lambda)',' T_{TE} / T_{TM} = 1');%,'location','northeast');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
set(leg,'box','on','linewidth',lwz)
tit=title('');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','HorizontalAlignment','right')
print('-depsc2',sprintf('TH_diam_a1_1.eps'), '-r300')

% Phase shift
% -----------
figure('name','null2')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   1400 600],'PaperUnits','points','PaperPosition',[400    0   2000 600]); % paper 
fig2 = subplot(1,4,1:2);
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[-.2 .1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('TE-TM phase mismatch','FontSize',fsz)
plot(lb_t,retard.*(2*pi)-pi,'k-o','linewidth',lwz)
plot(lb_t,lb_t.*0,'k--','linewidth',lwz)
leg=legend(' phase shift error  \epsilon(\lambda) (rad)',' \epsilon = 0');%,'location','northeast');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','southwest')
set(leg,'box','on','linewidth',lwz)
tit=title('N-band AGPM  [11-13.2 µm]  with \Lambda = 4.6 µm,  ');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')
% null depth
% -----------
fig5 = subplot(1,4,3:4);
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
%set(gca,'ylim',[170 195])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
ylabel('Null depth','FontSize',fsz)
plot(lb_t,nul_res_sp_b,'o-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,lb_t.*0+nuldpt,'k--','linewidth',2*lwz)
leg=legend(' Null depth N(\lambda)',' µ = mean null depth');%,'location','northeast');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
set(leg,'box','on','linewidth',lwz)
tit=title('F = 40%,  h = 13.5 µm  and  \alpha = 2.6°');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold','HorizontalAlignment','right')
print('-depsc2',sprintf('TH_diam_a1_2.eps'), '-r300')






break




















% % Calcul rcwa
% % -----------
% 
% vLb=[.45 .54 .99 1.3 1.68];
% vd=[4.94 5.94 11.19 14.79 19.77];
% vlb_min=[.53 .65 1.15 1.5 2];
% vlb_max=[.65 .79 1.4 1.8 2.4];
% lb0=vlb_min+(vlb_max-vlb_min)./2;
% F=0.61;%
% nlb=41;
% for i=1:5
%     Lb=vLb(i);
%     d=vd(i);
%     lb_min=vlb_min(i);
%     lb_max=vlb_max(i);
%     [fval1,fval2]=null
%     vlb(i,:)=lb_t;%./lb0(i);
%     nuldpt(i)=fval1;
%     sigerr(i)=fval2;
%     vnull(i,:)=nul_res_sp_b;
%     vretard(i,:)=retard;
% end 
% %x0=[d];
% %[x,fval,exitflag] = fminsearchbnd(@null,x0,d*.5,d*2,optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
% %x
% %[nuldpt,sigerr]=null(x)
% 
% % Figures
% % -------
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width 
% % All Phase shifts
% % ----------------
% close all
% i=4;
% figure('name','allps')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',[vlb_min(i) vlb_max(i)])%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[170 195])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
% ylabel('Phase shift  \Delta\Phi_{TE-TM} (°)','FontSize',fsz)
% plot(vlb(i,:),polyval(polyfit(vlb(i,:),vretard(i,:),7),vlb(i,:)).*360,'k-o','linewidth',lwz)
% % plot(vlb(2,:),polyval(polyfit(vlb(2,:),vretard(2,:),7),vlb(2,:)).*360,'k-o','linewidth',lwz)
% % plot(vlb(3,:),polyval(polyfit(vlb(3,:),vretard(3,:),7),vlb(3,:)).*360,'k-o','linewidth',lwz)
% % plot(vlb(4,:),polyval(polyfit(vlb(4,:),vretard(4,:),7),vlb(4,:)).*360,'k-o','linewidth',lwz)
% % plot(vlb(5,:),polyval(polyfit(vlb(5,:),vretard(5,:),7),vlb(5,:)).*360,'k-o','linewidth',lwz)
% plot(vlb(i,:),vlb(i,:).*0+180,'k--','linewidth',lwz)
% leg=legend(' \lambda = 1.5-1.8 µm (H band)',' \pi-phase shift');%,'location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% %print('-depsc2',sprintf('TH_sio2_ps4.eps'), '-r300');
% %' \lambda = 0.53-0.65 µm (visible)'
% %' \lambda = 0.65-0.78 µm (visible)'
% %' \lambda = 1.15-1.4 µm (J band)'
% %' \lambda = 1.5-1.8 µm (H band)'
% %' \lambda = 2-2.4 µm (K band)'
% 
% 
% break
% break
% 
% 
% % Calcul rcwa
% % -----------
% 
% Lb=1.68;%
% d=19.76;%
% F=0.61;%
% nlb=41;
% [nuldpt,sigerr]=null()
% 
% l=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
% TH_sio2 % Data absorption
% close all
% % ajustement courbe
% tmp1(1)=tmp1(1)*1.085;
% tmp1(2)=tmp1(2)*.93;
% tmp2(1)=tmp2(1)*1.03;
% tmp2(2)=tmp2(2)*.98;
% tmp2(3)=tmp2(3)*1.04;
% tmp2(4)=tmp2(4)*1.10;
% tmp2(5)=tmp2(5)*1.06;
% tmp2=polyval(polyfit(lb_t,tmp2,11),lb_t);
% tmp4(2)=tmp4(2)*.7;
% tmp5(2)=tmp5(2)*.8;
% tmp5(4)=tmp5(4)*1.5;
% tmp5(6)=tmp5(6)*.8;
% tmp5(15)=tmp5(15)*.5;
% tmp5(16)=tmp5(16)*.4;
% tmp5(21)=tmp5(21)*1.4;
% tmp1=pchip(lb_t,tmp1,lb_t);
% tmp2=pchip(lb_t,tmp2,lb_t);
% tmp4=pchip(lb_t,tmp4,lb_t);
% tmp5=pchip(lb_t,tmp5,lb_t);
% Tin =pchip(lb_t,(tmp1+tmp2)./2,lb_t);
% Rin =pchip(lb_t,(tmp4+tmp5)./2,lb_t);
% TNUM = T.*Tin.*abs;
% TDEN = 1-R.*Rin.*abs.^2;
% Ttot = TNUM ./ TDEN;
% RNUM = R.*Tin.^2.*abs.^2;
% RDEN = 1-R.*Rin.*abs.^2;
% Rtot = Rin + RNUM ./ RDEN;
% mean(Ttot)
% mean(Rtot)
% 
% % Figures
% % -------
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width 
% 
% 
% % Total Trans
% % -----------
% figure('name','trans')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[.82 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
% ylabel('Transmittance','FontSize',fsz)
% plot(lb_t,pchip(lb_t,T,lb_t),'--','color',myblue,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
% plot(lb_t,pchip(lb_t,Ttot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
% leg=legend(' T_{out} = (n-1)^2 / (n+1)^2',' T_{in} = (T_{TE}+T_{TM}) / 2',' T_{TM}',' T_{TE}',' Total transmittance');%,'location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','southeast')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('TH_sio2_Ttot.eps'), '-r300');
% 
% % Total Refl
% % ----------
% figure('name','relf')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[0 .085])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
% ylabel('Reflectance','FontSize',fsz)
% plot(lb_t,pchip(lb_t,Rtot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
% plot(lb_t,pchip(lb_t,R,lb_t),'--','color',myblue,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
% leg=legend(' Total reflectance',' R_{out} = (n-1)^2 / (n+1)^2',' R_{in} = (R_{TE}+R_{TM}) / 2',' R_{TM}',' R_{TE}');%,'location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('TH_sio2_Rtot.eps'), '-r300');

% 
% % Phase shift
% % -----------
% figure('name','ps')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[170 195])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
% ylabel('Phase shift  \Delta\Phi_{TE-TM} (°)','FontSize',fsz)
% %plot(lb_t,retard.*360,'k-o','linewidth',lwz)
% plot(lb_t,polyval(polyfit(lb_t,retard,7),lb_t).*360,'k-o','linewidth',lwz)
% plot(lb_t,lb_t.*0+180,'k--','linewidth',lwz)
% leg=legend(' HWP phase shift',' \pi-phase shift');%,'location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('TH_sio2_ps.eps'), '-r300');
% 
% 
% break 
% break
% 
% 
% % Optim dn form
% % -------------
% Lb=1.68;
% F=0.61;
% d=19.78;
% x0=d;
% nlb=13;
% [x,fval,exitflag] = fminsearchbnd(@null,x0,15,25,optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
% x
% F_min = 0.57;
% F_max = 0.65;
% nbF=5;
% for i=1:nbF
%     F=F_min+(F_max-F_min)*(i-1)/(nbF-1);
%     [nulldpt(i),sigerr(i)]=null(x);
%     retardtot(i,:)=retard;
% end
% % Figures
% % -------
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width
% close all
% figure('name','refl')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[.045 .075])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz)
% ylabel('Form birefringence   \Deltan_{TE-TM}','FontSize',fsz)
% plot(lb_t,polyval(polyfit(lb_t,retardtot(1,:),2),lb_t).*lb_t./d,':d','color',mygreen,'linewidth',2*lwz,'markersize',8)
% plot(lb_t,polyval(polyfit(lb_t,retardtot(2,:),2),lb_t).*lb_t./d,':+','color',mygreen,'linewidth',2*lwz,'markersize',8)
% plot(lb_t,polyval(polyfit(lb_t,retardtot(3,:),2),lb_t).*lb_t./d,':o','color',mygreen,'linewidth',2*lwz,'markersize',8)
% plot(lb_t,polyval(polyfit(lb_t,retardtot(4,:),2),lb_t).*lb_t./d,':x','color',mygreen,'linewidth',2*lwz,'markersize',8)
% plot(lb_t,polyval(polyfit(lb_t,retardtot(5,:),2),lb_t).*lb_t./d,':s','color',mygreen,'linewidth',2*lwz,'markersize',8)
% plot(lb_t,lb_t./2/d,'k','linewidth',lwz)
% leg=legend(' F = 0.57',' F = 0.59',' F = 0.61',' F = 0.63',' F = 0.65',' ideal')%,'location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',.9*fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('TH_sio2_dn.eps'), '-r300');
% 
% 
% 
% break
% break
% 
% 
% 
% % Optimisation (Lb,d)
% % -------------------
% disp('start')
% datestr(now)
% tic
% lb_min=.53;%.65;%1.15;%1.5;%2;%
% lb_max=.65;%.79;%1.4;%1.8;%2.2;%
% F=.6;
% d=5.3;%6.5;%11.5;%15;%20;%
% d_min=d/2;
% d_max=d*1.5;
% x0=[F d];
% Lb_min=.4;%.49;%0.8625;%1.125;%1.5;%
% Lb_max=.49;%.59;%1.05;%1.35;%1.8;%
% nbLb=11;%21;%5;%
% for i=1:nbLb
%     i
% %    F=F_min+(F_max-F_min)*(i-1)/(nbF-1);
%     Lb=Lb_min+(Lb_max-Lb_min)*(i-1)/(nbLb-1);
% %    [x,fval,exitflag] = fminsearchbnd(@null,x0,[.2 10],[.8 25],optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
%     [x,fval,exitflag] = fminsearchbnd(@null,x0,[.2 d_min],[.8 d_max],optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%     x
%     [fval1,fval2]=null(x)
%     Xaxis(i)=Lb;
%     Yaxis1(i)=x(1);
%     Yaxis2(i)=x(2);
%     nuldpt(i)=fval1;
%     sigerr(i)=fval2;
% end
% disp('end')
% datestr(now)
% toc
% % Figures
% % -------
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width
% band=[Lb_min Lb_max];
% % overplot
% % --------
% figure('name','refl')
% set(gcf,'color',[1 1 1])
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[1 1.6])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% %[AX,H1,H2] = plotyy(Xaxis,Yaxis1,Xaxis,Yaxis2);
% [AX,H1,H2] = plotyy(Xaxis,polyval(polyfit(Xaxis,Yaxis1,4),Xaxis)-.03,Xaxis,polyval(polyfit(Xaxis,Yaxis2,4),Xaxis)-1.1);
% %retmoy = zeros(size(lb_t)) + mean(retard);
% %plot(lb_t,retmoy,'k-.','Linewidth',2)
% %plot(lb_t,zeros(size(lb_t)) + .5,'k--')
% %axis tight
% set(AX(1),'FontSize',fsz,'FontWeight',fwz,'XLim',band,'XMinorTick','on','YColor','k')%,'ylim',[.56 .72])
% set(AX(2),'FontSize',fsz,'FontWeight',fwz,'XLim',band,'XMinorTick','on','YColor','k')%,'ylim',[17 21])
% %set(AX(1),'ylim',[.56 .72],'ytick',[.56 .58 .6 .62 .64 .66 .68 .7 .72])
% %set(AX(2),'ylim',[17 21],'ytick',[17 17.5 18 18.5 19 19.5 20 20.5 21])
% set(H1,'color',myblue,'marker','o','linewidth',lwz)%,'marker','.')
% set(H2,'color',mygreen,'marker','v','linewidth',lwz)%,'marker','.')
% %set([H1;H2],'Color','k','marker','o')
% xlabel('Period  \Lambda (µm)','FontSize',fsz)
% ylabel(AX(1),'Filling factor  F','FontSize',fsz)
% ylabel(AX(2),'Grating depth  h (µm)','FontSize',fsz)
% leg=legend(' Filling factor',' Depth','location','northeast');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',.8*fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% %print('-depsc2',sprintf('TH_sio2_over.eps'), '-r300');
% % nuldpt
% % ------
% figure('name','refl')
% set(gcf,'color',[1 1 1])
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',band)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[1 1.6])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% %plot(Xaxis,polyval(polyfit(Xaxis,sigerr,4),Xaxis),'k','linewidth',2,'marker','o')
% plot(Xaxis,nuldpt,'k','linewidth',2,'marker','o')
% xlabel('Period  \Lambda (µm)','FontSize',fsz)
% ylabel('Null depth','FontSize',fsz)
% Sigerr
% ------
% figure('name','refl')
% set(gcf,'color',[1 1 1])
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',band)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% %set(gca,'ylim',[1 1.6])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% plot(Xaxis,polyval(polyfit(Xaxis,sigerr,6),Xaxis),'k','linewidth',2,'marker','o')
% %plot(Xaxis,sigerr,'k','linewidth',2,'marker','o')
% xlabel('Period  \Lambda (µm)','FontSize',fsz)
% ylabel('Phase shift error std  \sigma','FontSize',fsz)
% % print('-depsc2',sprintf('TH_sio2_sigerr.eps'), '-r300');


% break
% save sio2_sigerr2.mat










