%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Libération
%----------
save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost bandAR PS_target lb_target ps_offset Retard_target

%Allocation mémoire
%------------------
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
p=1;%13;%3;%
%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=16;%25;%50;%8;%


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
EI_choice=1;%3;%
%Permittivités milieu extérieur émergent
EIII_choice=14;%1;%
%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=14;
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


% AFTA retarder
filter = 2;%3;%1;%
lb_vec = [.430 .550 .980];
Lb_vec = [.262 .338 .607];
lb = lb_vec(filter);
Lb = Lb_vec(filter);
ps_offset = 0;
PS_target_vec = mod([0.0356 0.0073 0.0136]*2*pi,2*pi)+ps_offset;
%PS_target_vec = pi+mod([0.0356 0.0073 0.0136]*2*pi,2*pi)+ps_offset;
PS_target_vec(PS_target_vec>pi) = PS_target_vec(PS_target_vec>pi)-2*pi;
PS_target = PS_target_vec(filter);



Dlbmm = [0.2 0.1 0.03];  %bandwidth
d=0.1414;%
F=0.8948;%
for mm = 1:3

% Band
% ------
nlb=41;%2;%
Dlb = Dlbmm(mm);%0.03;%0.2;
lb_min = lb*(1-Dlb/2);
lb_max = lb*(1+Dlb/2);
bandx = [lb_min lb_max];
bandAR = [lb_min (lb_max+lb_min)/2 lb_max];

% Absorption
% ----------
lb_t=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
TH_sio2_new

%Paramètres du réseau
%--------------------
% période --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
%d=0.146;%
%F=0.885;%
%Lb=.338;%.25;%.263;%.278;%.641;%.607;%
pente_deg=2;%
pente=deg2rad(pente_deg);
%COND LIMIT
Fmax = 1-2*d*tan(pente)/Lb
dmax = (1-F)*Lb/(2*tan(pente))

% Calcul rcwa
% -----------
%x0=[F d Lb];
if mm == 1
    x0=[d F];
    [x,fval,exitflag] = fminsearchbnd(@null_PS,x0,[0 0],[1 1],optimset('MaxIter',35,'Display','iter','TolX',1e-10,'TolFun',1e-10));
    d=x(1);
    F=x(2);
    %Lb=x(3);



% Figures
% -------
close all

% % PS error
% % --------
% newFig
% set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'ylim',[-300 600])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
% ylabel('Retardation (nm)','FontSize',fsz)
% %plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
% %plot(lb_t,lb_t*0+PS_target,'--','color',myblue,'linewidth',lwz)
% %plot(lb_t,abs(tmp3-PS_target)./PS_target.*100,'-','color',mygreen,'linewidth',2*lwz)
% plot(lb_t,tmp3./(2*pi).*lb_t.*1e3,'-','color',mygreen,'linewidth',2*lwz)
% %plot(lb_t,lb_t*0+Retard_target,'--','color',myblue,'linewidth',lwz)
% plot(lb_t,lb_t./(2*pi).*PS_target.*1e3,'--','color',myblue,'linewidth',lwz)
% title(sprintf('Period = %3.3g nm, F = %.2g, h = %d nm',Lb*1e3,F,d*1e3))
% leg=legend(' retarder',' target');
% set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
% tick2latex
% 
% %save
% print('-depsc2',sprintf('sio2_retard_lb%d_Dlb=%.2g.eps',filter,Dlb), '-r300');


% break

% Phase TE-TM
% -----------
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[-2 0.5])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Waves','FontSize',fsz)
%plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
plot(lb_t,tmp3./(2*pi),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,lb_t*0+PS_target./(2*pi),'--','color',myblue,'linewidth',lwz)
%plot(lb_t,tmp3-PS_target,'-','color',mygreen,'linewidth',2*lwz)
title(sprintf('Period = %3.3g nm, F = %.2g, h = %d nm',Lb*1e3,F,d*1e3))
leg=legend(' retarder',' target');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex

%save
print('-depsc2',sprintf('sio2_ps_lb%d_Dlb=%.2g.eps',filter,Dlb), '-r300');


else
    x0=[d];
    %[x,fval,exitflag] = fminsearch(@null_PS_optim,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-10,'TolFun',1e-10));
    [x,fval,exitflag] = fminsearchbnd(@null_PS_optim,x0,[0],[1],optimset('MaxIter',15,'Display','iter','TolX',1e-10,'TolFun',1e-10));
    d=x(1);
end
F
dmm(mm)=d
psmm(mm)= fval
ps2mm(mm)= fval/(2*pi)  % waves
end


break


% % new routine
% delta_lb=0.1; %band
% nlb = 11;%
% nbF = 8;
% F_min = 0.1;
% F_max = 0.9;
% lb_target_temp = lb_target;
% PS_target_temp = PS_target;
% d=10;
% x0=[d];
% for kk=3
%     lb_min=(1-delta_lb/2)*lb_target(kk);
%     lb_max=(1+delta_lb/2)*lb_target(kk);
%     lb_target(1:nlb)=lb_target_temp(kk);
%     PS_target(1:nlb)=PS_target_temp(kk);
%     for ii=1:nbF
%         kk
%         ii
%         F = F_min+(F_max-F_min)/(nbF-1)*(ii-1);
% %        [x,fval,exitflag] = fminsearchbnd(@null_PS,x0,Lb/5,Lb*10,optimset('MaxIter',35,'Display','iter','TolX',1e-10,'TolFun',1e-10));
%         [x,fval,exitflag] = fminsearch(@null_PS,x0,optimset('MaxIter',35,'Display','iter','TolX',1e-10,'TolFun',1e-10));
%         xParam(kk,ii) = F;
%         yParam(kk,ii) = x;
%     end
% end
% 
% newFig
% %set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% xlabel('Fill factor','FontSize',fsz)
% ylabel('Depth','FontSize',fsz)
% plot(xParam(1,:),yParam(1,:),'-','color',mygreen,'linewidth',2*lwz)
% plot(xParam(2,:),yParam(2,:),'-','color',myblue,'linewidth',2*lwz)
% plot(xParam(3,:),yParam(3,:),'-','color',myred,'linewidth',2*lwz)
% 
% break


% Null Depth + GHOST
% ------------------
close all
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Null depth','FontSize',fsz)
plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
plot(lb_t,nul_res_sp_b+NARGghost,'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,nul_res_sp_b,'--','color',myblue,'linewidth',1.2*lwz)
%str = sprintf('Outer L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('\Lambda = %3.3g \mu m, h = %3.3g \mu m,  F = %.2g',Lb,d,F);
%str = sprintf('Central L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%title(str,'FontSize',fsz)
leg=legend(' without ARG',' with ARG',' no ghost');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
tick2latex


%save
%print('-depsc2',sprintf('agpm_L_outer.eps'), '-r300');
print('-depsc2',sprintf('agpm_L_Null.eps'), '-r300');
%print('-dpng',sprintf('agpm_L_brun.png'), '-r300');





% Total Trans WITH ARG
% --------------------
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'ylim',[.75 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
plot(lb_t,pchip(lb_t,TARG,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,TARGtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
%tit=title('L-band AGPM with ARG');
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{out}} = T_{\rm{ARG}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','south')
tick2latex
print('-depsc2',sprintf('agpm_L_TtotARG.eps'), '-r300');










