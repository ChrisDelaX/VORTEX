%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Liberation
%----------
save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret 
global d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice 
global E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice 
global theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost 
global bandAR bandTAR nul_ghost nul_ghost_ARG nul_ideal
global tmp1p1 tmp1m1 tmp2p1 tmp2m1

%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)
N=8;%12;%
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
%18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, 
%    sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trapézoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trapézoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, 
%    sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)
p=13;%1;%3;%
%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=25;%50;%8;%


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
E_man=1.494;%1.44;
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


%Onde incidente
%--------------
%Angle incidence non conique
theta=deg2rad(0); theta_min=deg2rad(30); theta_max=deg2rad(50);
%Angle incidence conique
phi=deg2rad(0); phi_min=deg2rad(0); phi_max=deg2rad(0);
%Angle polarisation
psi=deg2rad(45); psi_min=deg2rad(0); psi_max=deg2rad(0);


% Parametres du reseau
%--------------------
% period --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391

nlb=41;%11;%
name = 'AGPM-METIS-L';
switch name
    case 'AGPM-METIS-L'
        lam_range = [2.9 4.1];
        AR_range = [3.2 3.8];
        Lb=1.21; %period
        lw = 0.5; % line width
        dlw = 0.05;
        dep = 4.4; % depth
        ddep = 0.05;
        sw = 3; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N-BT2-new3'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.7; % line width
        dlw = 0.05;
        dep = 14.5; % depth
        ddep = 0.;
        sw = 3.1; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N-BT2-new2'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.7; % line width
        dlw = 0.05;
        dep = 15.5; % depth
        ddep = 0.15;
        sw = 3.1; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N-BT2'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.7; % line width
        dlw = 0.05;
        dep = 15.5; % depth
        ddep = 0.15;
        sw = 3.5; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N-BT1-new3'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.8; % line width
        dlw = 0.04;
        dep = 13.9; % depth
        ddep = 0.0;
        sw = 2.6; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N-BT1-new2'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.8; % line width
        dlw = 0.0;
        dep = 14.; % depth
        ddep = 0.15;
        sw = 2.6; % sidewall angle
        dsw = 0.1;        
    case 'AGPM-N-BT1'
        lam_range = [10 12.5];
        AR_range = lam_range;
        Lb=4; %period
        lw = 1.8; % line width
        dlw = 0.05;
        dep = 14.9; % depth
        ddep = 0.15;
        sw = 3; % sidewall angle
        dsw = 0.1;
    case 'AGPM-N3'
        lam_range = [10 12.5];
        AR_range = [12.5 13.3];
        Lb=4.6; %period
        lw = 1.66; % line width
        dlw = 0.05;
        dep = 12.6; % depth
        ddep = 0.15;
        sw = 2.75; % sidewall angle
        dsw = 0.1;
    case 'opti_F.45_N'
        lam_range = [11 13.2];
        AR_range = lam_range;
        Lb=4.6; %period
        lw = Lb*.45; % line width
        dlw = 0.05;
        dep = 16.9; % depth
        ddep = 0.15;
        sw = 3; % sidewall angle
        dsw = 0.1;
    case 'opti_F.45_L'
        lam_range = [3.5 4.1];
        AR_range = lam_range;
        Lb=1.42; %period
        lw = Lb*.45; % line width
        dlw = 0.05;
        dep = 5.21; % depth
        ddep = 0.15;
        sw = 3; % sidewall angle
        dsw = 0.1;        
end
x0=[];

% Absorption
% ----------
lb_min = lam_range(1);
lb_max = lam_range(2);
lb_t=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
bandx = [lb_min lb_max];
TH_diam_new


% Anti-reflection grating
% -----------------------
bandAR = [AR_range(1) (AR_range(2)+AR_range(1))/2 AR_range(2)];
bandTAR = [.984 .996 .984];


% Calculate Null Depth
% --------------------
% on spec: sidewall ON
pente1 = sw
pente = deg2rad(pente1);
d=dep;
F=lw/Lb;
[fval1,fval2] = null(x0);
nul_on1 = nul_res_sp_b;
% off max
d=dep+ddep;
F=(lw+dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmax1 = nul_res_sp_b;
% off min
d=dep-ddep;
F=(lw-dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmin1 = nul_res_sp_b;

% on spec: sidewall MAX
pente2 = sw+dsw
pente = deg2rad(pente2);
d=dep;
F=lw/Lb;
[fval1,fval2] = null(x0);
nul_on2 = nul_res_sp_b;
% off max
d=dep+ddep;
F=(lw+dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmax2 = nul_res_sp_b;
% off min
d=dep-ddep;
F=(lw-dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmin2 = nul_res_sp_b;

% on spec: sidewall MIN
pente3 = sw-dsw
pente = deg2rad(pente3);
d=dep;
F=lw/Lb;
[fval1,fval2] = null(x0);
nul_on3 = nul_res_sp_b;
% off max
d=dep+ddep;
F=(lw+dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmax3 = nul_res_sp_b;
% off min
d=dep-ddep;
F=(lw-dlw)/Lb;
[fval1,fval2] = null(x0);
nul_offmin3 = nul_res_sp_b;

%---------

close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width


band = lam_range;
%bandtick = [10.0 10.5 11.0 11.5 12.0 12.5];
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];


figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); 
% Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',band)%,'xtick',bandtick)
set(gca,'YScale','log')
set(gca,'YLim',[1e-4 1])
%set(gca,'ylim',[.0001 .2])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('Null Depth','FontSize',fsz*1.2)
plot(lb_t,nul_offmax2,'g--','linewidth',lwz,'HandleVisibility','off');
y1=plot(lb_t,nul_on2,'g','linewidth',lwz);
plot(lb_t,nul_offmin2,'g:','linewidth',lwz,'HandleVisibility','off');
plot(lb_t,nul_offmax1,'k--','linewidth',lwz,'HandleVisibility','off');
y1=plot(lb_t,nul_on1,'k','linewidth',lwz);
plot(lb_t,nul_offmin1,'k:','linewidth',lwz,'HandleVisibility','off');
plot(lb_t,nul_offmax3,'r--','linewidth',lwz,'HandleVisibility','off');
y1=plot(lb_t,nul_on3,'r','linewidth',lwz);
plot(lb_t,nul_offmin3,'r:','linewidth',lwz,'HandleVisibility','off');

l=legend(sprintf(' wall angle = %3.2f°',pente2),  sprintf(' wall angle = %3.2f°', ...
    pente1), sprintf(' wall angle = %3.2f°',pente3));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)

t=title(sprintf('p=%3.2fµm, d=%3.2f±%3.2fµm, lw=%3.2f±%3.2fµm, AR:[%3.2f-%3.2fµm]', ...
    Lb,dep,ddep,lw,dlw,AR_range(1),AR_range(2)));
set(t,'Fontname',fnz,'FontSize',fsz*.8,'FontWeight',fwz,'HorizontalAlignment','center')

print('-dpng', [name,'.png'], '-r300');



















