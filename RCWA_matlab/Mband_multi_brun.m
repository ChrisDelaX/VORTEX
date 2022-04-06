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
Lb=0.8;%1.42;%
Lb_min=4.7;Lb_max=4.75;

%Epaisseurs µm
d=4.2;%4.3;%
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%Facteurs de remplissage
F=.45;%0.4;%

%Onde incidente
%--------------

%Longueur d'onde µm
lb=2.15;%11;%      35 9.5 738
lb_min=4.6;%1.95;%1.48;%
lb_max=5.0;%2.35;%
nlb=11;%81;%31;%

%Angle incidence non conique
theta=rad(0);
theta_min=rad(30); theta_max=rad(50);
%Angle incidence conique
phi=rad(0);
phi_min=rad(0); phi_max=rad(0);
%Angle polarisation
psi=rad(45);
psi_min=rad(0); psi_max=rad(0);


% ABSORPTION
% ----------

l=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
TH_diam 
close all

% Calcul rcwa
% -----------

pente_min = 2.5;%2.75;%
pente_max = 3.5;%3.25;%
Lb = 1.89;
F_min=0.30;%0.4;%0.36;%
F_max=0.50;%0.5;%0.48;%
d_min=4.5;%2.5;%2.4;%
d_max=7.5;%3.6;%3.8;%

npts_pente = 1;%3;%5;%11;%
npts_F = 24;%4;%8;%28;%42;%2;%1440;%
npts_d = 24;%4;%8;%18;%27;%2;%878;%

%spécial: Lb vs F
%d=5.3;%4.8;%5.9;%5.2;%4.5;%
pente_deg = 3;%3.2;
pente=rad(pente_deg);
Lb_min=2*d*tan(pente)/(1-F_max);%-.01;
Lb_max = 1.41;%;%1.5;
Lb_min = 1.09;
npts_Lb = npts_d;

tic
depart_tmp=now;
depart=datestr(depart_tmp)
for i=1:npts_pente
    %pente_deg = pente_min+(pente_max-pente_min)*(i-1)/(npts_pente-1)
    pente=rad(pente_deg);
    for j=1:npts_F
        j
        F = F_min+(F_max-F_min)*(j-1)/(npts_F-1);
        Xparam(j)=F;
        for k=1:npts_d
            d = d_min+(d_max-d_min)*(k-1)/(npts_d-1);
            Yparam(k)=d;
            %Lb = Lb_min+(Lb_max-Lb_min)*(k-1)/(npts_Lb-1);
            %Xparam(k)=Lb;
            [fval1,fval2] = null();
            %fval1 = 0.001*max(j,k);
            if imag(fval1) == 0
                nuldpt(j,k)=fval1;
            else
                nuldpt(j,k)=1;
            end
            Reject(j,k)=1/nuldpt(j,k);
            if i==1 && j==1 && k==1
                fin=datestr(depart_tmp+toc/3600/24*npts_pente*npts_F*npts_d)
            end
        end
    end
    fmulti = nuldpt;    
    dessin_multi
end
depart
fin
toc
datestr(now)
save('backup')

%load Lb.vs.F_depth=5.60_angle=3.00.mat

% Figures
% -------
close all
figure
hold on
%grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.4,0.4,0.4];
fnz = 'Arial'; % fontname: Helvetica
fsz = 16; % fontsize: 10
fwz = 'normal';  % fontweight: bold
msz = 8; % marker size
lwz = 2;  % line width
set(gcf,'color',mywhite)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'ticklength',-.9*get(gca,'ticklength'))
%set(0,'DefaultTextInterpreter', 'latex')
set(gcf(), 'Renderer', 'painters')

Zparam = log10(fmulti');
colorbarLim = [min(min(Zparam)) -0.9];
Zparam(Zparam>colorbarLim(2)) = colorbarLim(2);
hS = surfc(Xparam,Yparam,Zparam);
shading('interp') %flat, faceted
set(hS,'EdgeColor', 'none')
hC = colorbar;
set(hC,'box','on','linewidth',lwz/2)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(hC,'YLim',colorbarLim)

%contours
[C01,C02] = contour(Xparam,Yparam,Zparam,[-3.8:.2:-3.0],'k-');
[C03,C04] = contour(Xparam,Yparam,Zparam,[-2.7:.5:-1.7],'k--');
%[C05,C06] = contour(Xparam,Yparam,Zparam,[-1.6:.4:-0],'k-.');
clabel(C01,C02,'FontSize',fsz);
clabel(C03,C04,'FontSize',fsz);
clabel(C05,C06,'FontSize',fsz);

%graphics
% [V_nul,I_nul] = min(nuldpt);
% k=0;
% while log10(V_nul(k+1))>=-3.2
%     k=k+1;
% end
% best_Lb_min=Yparam(k);
% best_Lb_max=1.42;%
% best_Lb_rat=best_Lb_min/best_Lb_max;
% [idx idx] = min(abs(Yparam-best_Lb_max)); %index of closest value
% best_F_min=Xparam(I_nul(idx)); %closest value 
% best_F_max=Xparam(I_nul(k));
% best_F_moy=(best_F_max+best_F_min)/2; 
% line([best_F_moy-.01 best_F_moy-.01],[best_Lb_max best_Lb_min],'color',mywhite,'LineStyle','--','linewidth',lwz)
% line([best_F_moy+.01 best_F_moy+.01],[best_Lb_max best_Lb_min],'color',mywhite,'LineStyle','--','linewidth',lwz)
% line([best_F_moy best_F_moy],[best_Lb_max Lb_min],'color',mywhite,'linewidth',lwz)
% line([best_F_moy best_F_moy+.01],[best_Lb_max best_Lb_max],'color',mywhite,'LineStyle','--','linewidth',lwz)
% line([best_F_moy best_F_moy+.01],[best_Lb_min best_Lb_min],'color',mywhite,'LineStyle','--','linewidth',lwz)
% line([F_min best_F_moy],[best_Lb_max best_Lb_max],'color',mywhite,'linewidth',lwz)
% line([F_min best_F_moy],[best_Lb_min best_Lb_min],'color',mywhite,'linewidth',lwz)
axis([Xparam(1) Xparam(end) Yparam(1) Yparam(end)])
%xlabel('Period $\Lambda (\mu m)$','FontSize',fsz)
%ylabel('Filling factor $F$','FontSize',fsz)
%set(get(hC,'YLabel'),'String','Null Depth (log$_{10}$)','Rotation',-90,'FontSize',fsz)
%title('L band $[3.5 - 4.1\mu m]$','FontSize',fsz*1.5)
xlabel('Filling factor','FontSize',fsz*1.2)
ylabel('Depth (µm)','FontSize',fsz*1.2)
title('M band [4.6 - 5.0 µm] - Null Depth (log_{10})','FontSize',fsz*1.5)
%set(gca,'XTick',[1.1 1.15 1.2 1.25 1.3 1.35 1.4])
%set(gca,'YTick',[.41 .43 .45 .47 .49 .51 .53])
%set(gca,'YTick',[ ])
set(hC,'YTick',[-3.2 -3 -2.8 -2.6 -2.4 -2.2 -2 -1.8 -1.6 -1.4 -1.2 -1])
%str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: \\alpha = %1.0g^{\\circ}$$.  For $$h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %3.2g \\rightarrow F = %3.2g $$',best_Lb_max,pente_deg,d,best_Lb_rat,best_F_moy);
%tit=title(str,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold');

%tick2latex
%colorbar2latex(hC)


%save
%save(sprintf('Lb.vs.F_depth=%3.2f_angle=%3.2f.mat',d,pente_deg))

%print('-depsc2',sprintf('Lb.vs.F_depth=%3.2f_angle=%3.2f.eps',d,pente_deg), '-r300')
print('-depsc2',sprintf('d.vs.F_period=%3.2f_angle=%3.2f.eps',Lb,pente_deg), '-r300')
print('-dpng',sprintf('d.vs.F_period=%3.2f_angle=%3.2f.png',Lb,pente_deg), '-r300')

break
close all


%Full screen
pos0=get(0,'Screensize');
pos1=get(gcf,'Position');
pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



