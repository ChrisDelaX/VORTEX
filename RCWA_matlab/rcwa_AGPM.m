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

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1p1 tmp1m1 tmp2p1 tmp2m1 tmp1 tmp2 tmp3 tmp4 tmp5
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
L=25;%100;%

pente_deg=2.75;%0;%
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

E_man=1.5; % si manuel (8)

%Permittivités milieu extérieur incident
EI_choice=1;
%Permittivités milieu extérieur émergent
EIII_choice=3;
%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;

%Permittivités de la couche antireflet
E_AR_choice=1;%14;%
E_AR1_choice=1;
E_AR2_choice=1;
%Permittivités de la couche d'arret
E_arret_choice=1;


%Paramètres du réseau
%--------------------

%Pas du réseau
Lb=4.6;
%Epaisseurs
d=13.69;
%Facteurs de remplissage
F=0.4;


%Onde incidente
%--------------

%Longueur d'onde
lb=11;
lb_min=11;
lb_max=13.2;
nlb=11;%81;%

%Angle incidence non conique
theta=rad(0);
%Angle incidence conique
phi=rad(0);
%Angle polarisation
psi=rad(45);


% Calcul Null / RMS_err  (optimistation)
%---------------------------------------

disp('start')
datestr(now)
tic

pente_min = 0;%2.5;%
pente_max = 2.75;%3;%
npts_pente = 11;%3;%

for k=1:npts_pente    
    k
    pente_deg = pente_min+(pente_max-pente_min)*(k-1)/(npts_pente-1);
    pente=rad(pente_deg);
    %F_equiv=0.4+d*tan(pente)/Lb

    [fval1,fval2] = null();  % calculate Null Depth and Phase Shift Error RMS
    nuldpt = fval1;
    rms_err = fval2;
    save(['DATA_slope=',sprintf('%3.2f',pente_deg),'_F=',sprintf('%3.2f',F),'.mat'])
end


disp('end')
datestr(now)
toc



% Représentation graphique
% ------------------------
% ------------------------

close all
warning off MATLAB:polyfit
fnz = 'Arial'; % fontname
fsz = 20; % fontsize
fwz = 'Normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width
band = [11 13.2];
bandtick = [11.0 11.5 12.0 12.5 13.0];
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];


figure('name','transmittance')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'ylim',[60 100],'ytick',[60 70 80 90 100])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on','YMinorGrid','on')
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('Transmittance (%)','FontSize',fsz*1.2)

for k=1:npts_pente
    pente_deg = pente_min+(pente_max-pente_min)*(k-1)/(npts_pente-1);
    pente=rad(pente_deg);
    load(['DATA_slope=',sprintf('%3.2f',pente_deg),'_F=',sprintf('%3.2f',F),'.mat'])
    Trans = 100.*(tmp1+tmp2)./2;
    
    % b=plot(lb_t,Trans,'color','k');   % non lissé
    b=plot(liss_lb_t,polyval(polyfit(lb_t,Trans,8),liss_lb_t),'color','k');   % lissé
    
    if pente==rad(0)
        set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
    elseif pente==rad(1.375)
        set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
    elseif pente==rad(2.75)
        set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    end
end


TsansZOG= n2lb.*(2./(n2lb+1)).^2;
plot(lb_t,100.*TsansZOG,'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[12 11 7 6 2 1];
idy=[1 2 3 4 5 6];
str={' alpha = 0°' ' ...' ' alpha = 1.375°' ' ...' ' alpha = 2.75°' ' without zog'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','southeast')
print('-dpng', ['transmittance.png'], '-r300');






