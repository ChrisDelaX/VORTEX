%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Gestion de la mémoire
%---------------------
%---------------------

%Libération
%----------

%save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof

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

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=16;%25;%50;%8;%
pente_deg=3.2;%2.75;% 
pente=rad(pente_deg);

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
EIII_choice=3;

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;

%Permittivités de la couche antireflet
E_AR_choice=1;%14;
E_AR1_choice=1;
E_AR2_choice=1;

%Permittivités de la couche d'arret
E_arret_choice=1;


%Paramètres du réseau
%--------------------

%Pas du réseau (sublambda!)
% --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
Lb=4.6;%1.42;%
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
lb=3.5;%11;%      35 9.5 738
lb_min=2;%1.5;%.65;%11;%3.25;%3.5;%
lb_max=2.4;%1.8;%.78;%13.2;%4.25;%4.1;%
nlb=81;%11;%31;%

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

%[nuldpt,sigerr]=null()


% % AGPM-N3
% % -------
% lb_min=11;%
% lb_max=13.2;%
% Lb=4.6;%
% nlb=41;
% F=0.3587;%0.44;%
% d=12.6;%5.8;%6;%6;%
% penN3=2.75;%3.25;%2.9;%pen_max;%
% pente=rad(penN3);
% % calcul réjection
% [fval1,fval2] = null();
% nul_tmp3 = nul_res_sp_b;
% FN3=F;
% hN3=d;
% % AGPM-N4
% % -------
% lb_min=11;%
% lb_max=13.2;%
% Lb=4.6;%
% nlb=41;
% F=0.3967;%0.41;%
% d=13.8;%4.7;%
% penN4=2.75;%3.0;%2.9;%
% pente=rad(penN4);
% % calcul réjection
% [fval1,fval2] = null();
% nul_tmp4 = nul_res_sp_b;
% FN4=F;
% hN4=d;
% % AGPM-L3
% % -------
% lb_min=3.5;%
% lb_max=4.1;%
% Lb=1.42;%
% nlb=41;
% F=0.44;%
% d=5.8;%6;%6;%
% penL3=3.25;%2.9;%pen_max;%
% pente=rad(penL3);
% % calcul réjection
% [fval1,fval2] = null();
% nul_tmp3 = nul_res_sp_b;
% FL3=F;
% hL3=d;
% AGPM-L4
% -------
lb_min=3.5;%
lb_max=4.1;%
Lb=1.42;%
nlb=41;
F=0.41;%
d=4.7;%
penL4=3.0;%2.9;%
pente=rad(penL4);
% calcul réjection
[fval1,fval2] = null();
nul_tmp4 = nul_res_sp_b;
FL4=F;
hL4=d;


% ABSORPTION
% ----------
l=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
TH_diam 
close all


% Figures
% -------
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
fnz = 'Arial'; % fontname
fsz = 26; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2.2;  % line width 
bandx = [lb_min lb_max];


% WITHOUT ARG
% -----------
%in:AGPM
Tin =pchip(lb_t,(tmp1+tmp2)./2,lb_t);
Rin =pchip(lb_t,(tmp4+tmp5)./2,lb_t);
%tot
TNUM = T.*Tin.*abs;
TDEN = 1-R.*Rin.*abs.^2;
Ttot = TNUM ./ TDEN;
RNUM = R.*Tin.^2.*abs.^2;
RDEN = 1-R.*Rin.*abs.^2;
Rtot = Rin + RNUM ./ RDEN;
TnoARG=mean(Ttot);
RnoARG=mean(Rtot);

% Total Trans
% -----------
close all
figure('name','trans')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
set(gca,'ylim',[.45 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,T,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Ttot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
leg=legend(' T_{TM}',' T_{TE}',' T_{in} = (T_{TE}+T_{TM}) / 2',' T_{out} = (n-1)^2 / (n+1)^2',' Total transmittance');%,'location','northeast');
tit=title('AGPM-L4 w/o ARG');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
print('-depsc2',sprintf('TH_diam_L4_Ttot.eps'), '-r300');


% Total Refl
% ----------
figure('name','relf')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Reflectance','FontSize',fsz)
set(gca,'ylim',[0 .4])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,Rtot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,R,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
leg=legend(' Total reflectance',' R_{out} = (n-1)^2 / (n+1)^2',' R_{in} = (R_{TE}+R_{TM}) / 2',' R_{TE}',' R_{TM}');%,'location','northeast');
tit=title('AGPM-L4 w/o ARG');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
print('-depsc2',sprintf('TH_diam_L4_Rtot.eps'), '-r300');


% WITH ARG
% --------

%out:ARG
ltmp=[3.5 3.8 4.1];
TARG=[.97 .99 .97];
%ltmp=[11 12.1 13.2];
%TARG=[.997 .998 .994];
TARG=polyval(polyfit(ltmp,TARG,2),l);
RARG=1-TARG;
%in:AGPM
Tin =pchip(lb_t,(tmp1+tmp2)./2,lb_t);
Rin =pchip(lb_t,(tmp4+tmp5)./2,lb_t);
%tot
TNUM = TARG.*Tin.*abs;
TDEN = 1-RARG.*Rin.*abs.^2;
TARGtot = TNUM ./ TDEN;
RNUM = RARG.*Tin.^2.*abs.^2;
RDEN = 1-RARG.*Rin.*abs.^2;
RARGtot = Rin + RNUM ./ RDEN;
TwithARG=mean(TARGtot);
RwithARG=mean(RARGtot);


% Total Trans
% -----------
close all
figure('name','trans')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
set(gca,'ylim',[.7 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,TARG,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,TARGtot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
leg=legend(' T_{out} = T_{ARG}',' T_{in} = (T_{TE}+T_{TM}) / 2',' T_{TM}',' T_{TE}',' Total transmittance');%,'location','northeast');
tit=title('AGPM-L4 with ARG');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
print('-depsc2',sprintf('TH_diam_L4_TtotARG.eps'), '-r300');


% Total Refl
% ----------
figure('name','relf')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Reflectance','FontSize',fsz)
set(gca,'ylim',[0 .2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,RARGtot,lb_t),'o-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,RARG,lb_t),'--','color',myblue,'linewidth',lwz)
leg=legend(' Total reflectance',' R_{TE}',' R_{TM}',' R_{in} = (R_{TE}+R_{TM}) / 2',' R_{out} = R_{ARG}');%,'location','northeast');
tit=title('AGPM-L4 with ARG');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
print('-depsc2',sprintf('TH_diam_L4_RtotARG.eps'), '-r300');



% GHOST
% -----
Tghost = Ttot-T.*Tin.*abs;
Nghost = Tghost./Ttot;
TARGghost = TARGtot-TARG.*Tin.*abs;
NARGghost = TARGghost./TARGtot;

figure('name','ghost')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength   \lambda  (\mum)','FontSize',fsz)
ylabel('Null depth','FontSize',fsz)
plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
plot(lb_t,nul_res_sp_b+NARGghost,'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,nul_res_sp_b,'-','color',myblue,'linewidth',lwz)
leg=legend(' Null + ghost w/o ARG',' Null + ghost with ARG',' Null (w/o ghost)');%,'location','northeast');
tit=title('AGPM-L4: Null depth + ghost');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
set(leg,'box','on','linewidth',lwz)
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
print('-depsc2',sprintf('TH_diam_L4_ghost.eps'), '-r300');

NullNO=mean(nul_res_sp_b);
NullYES=mean(nul_res_sp_b+Nghost);
NullYESARG=mean(nul_res_sp_b+NARGghost);

Tghostmean=mean(Tghost);
TARGghostmean=mean(TARGghost);
Nghostmean=mean(Nghost);
NARGghostmean=mean(NARGghost)

break 
break







