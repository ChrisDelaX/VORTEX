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

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost


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

% Kband
% -------
lb_min=1.95;%3.5;%3.5;%1.48;%
lb_max=2.35;%4.1;%4.1;%1.77;%
%lb_min2=1.95;%
%lb_max2=2.35;%
nlb=11;%81;%


% ABSORPTION
% ----------
l=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
TH_diam
bandx = [lb_min lb_max];
close all


Lb=1.42;%0.61;%0.8;%
F=0.4506;%0.45;%
d=5.2227;%1.9;%2.22;%2.5;%2.9368;%
x0=[F d];
penK=3.0;%2.5;%3.25;%2.75;%
pente=rad(penK);
% calcul réjection
[fval1,fval2] = null();
%[x,fval,exitflag] = fminsearchbnd(@null,x0,d*.5,d*2,optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%x




% Figures
% -------
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.4,0.4,0.4];
fnz = 'Arial'; % fontname: Helvetica
fsz = 12; % fontsize: 10
fwz = 'normal';  % fontweight: bold
msz = 8; % marker size
lwz = 2;  % line width
set(gcf,'color',mywhite)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(0,'DefaultTextInterpreter', 'latex')



% WITHOUT ARG
% -----------

% Total Trans
% -----------
close all
figure('name','trans')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'ylim',[.5 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,T0,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Ttot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
tit=title('H-band AGPM w/o ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{out}} = \frac{4n}{(n+1)^2}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_K_Ttot.eps'), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



% Total Refl
% ----------
figure('name','relf')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'ylim',[0 .5])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Reflectance','FontSize',fsz)
plot(lb_t,pchip(lb_t,Rtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,R0,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
tit=title('H-band AGPM w/o ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' Total refl.',' $R_{\rm{out}} = \frac{(n-1)^2}{(n+1)^2}$',' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_K_Rtot.eps'), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]




% WITH ARG
% --------

% Total Trans
% -----------
figure('name','trans')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'ylim',[.75 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
plot(lb_t,pchip(lb_t,TARG,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,TARGtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
tit=title('H-band AGPM with ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{out}} = T_{\rm{ARG}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_K_TtotARG.eps'), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



% Total Refl
% ----------
figure('name','relf')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Reflectance','FontSize',fsz)
set(gca,'ylim',[0 .25])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,RARGtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,RARG,lb_t),'--','color',myblue,'linewidth',lwz)
tit=title('H-band AGPM with ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' Total refl.',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$',' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',' $R_{\rm{out}} = R_{\rm{ARG}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_K_RtotARG.eps'), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



% GHOST
% -----

figure('name','ghost')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',lwz)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[1e-4 1e-2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Null depth','FontSize',fsz)
plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
plot(lb_t,nul_res_sp_b+NARGghost,'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,nul_res_sp_b,'-','color',myblue,'linewidth',lwz)
tit=title('H-band AGPM: Null depth + ghost');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' Null + ghost w/o ARG',' Null + ghost with ARG',' Null (w/o ghost)');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_optimal.eps'), '-r300');

%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]




NullNO=mean(nul_res_sp_b)
NullYES=mean(nul_res_sp_b+Nghost);
NullYESARG=mean(nul_res_sp_b+NARGghost)


Nghostmean=mean(Nghost)
NARGghostmean=mean(NARGghost);







