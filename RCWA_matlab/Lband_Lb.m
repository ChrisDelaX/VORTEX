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
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost bandAR


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
p=13;%1;%3;%
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
EIII_choice=3;%1;%
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


% Band
% ------
lb_min=3.5;%1.5;%2;%3.5;%4.6;%8;%11;%
lb_max=4.1;%2.4;%13.5;%5.0;%11.3;%13.2;%
%lb_min2=1.95;%
%lb_max2=2.35;%
lb_step = 0.025;
nlb=25;%251;%101;%17;%(lb_max-lb_min)/lb_step+1;%31;%46;%81;%11;%
%bandx = [3.5 4.1];
bandx = [lb_min lb_max];
ARcenter = 3.8;%1.9;%2.2;%4.1;
ARwidth = ARcenter/5;%/6;%%0.6;%0.7;%
bandAR = [ARcenter-ARwidth/2 ARcenter ARcenter+ARwidth/2];
%bandAR = [lb_min (lb_max+lb_min)/2 lb_max];


% Absorption
% ----------
lb_t = lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
%lb_t = lb_t(lb_t<=4.1 | lb_t>=4.6);
%lb_t = lb_t(lb_t<=1.8 | lb_t>=2.0);
TH_diam_new


%Paramètres du réseau
%--------------------
% période --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
Lb=1.42;%0.6893;%1.4;%1.5349;%1.5;%1.89;%3.32;%4.58;%
F=0.45;%0.66;%0.45;%0.4766;%0.47;%.44;%    
d=5.2;%3.3393;%5.5;%7.0638;%7;%5;%  
%x0=[d];
x0=[F d Lb];
%x0=[F d];
pente_deg=2.8;%3.15;%.15;%3.2;%2.5;%3.25;%2.75;%
pente=deg2rad(pente_deg);
        
%COND LIMIT
%Fmax = 1-2*d*tan(pente)/Lb
%dmax = (1-F)*Lb/(2*tan(pente))
%break
%
% Calcul rcwa
% -----------
%[fval1,fval2] = null();
[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2 Lb/2],[F*2 d*2 Lb*2],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,F/2,F*2,optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,d/2,d*2,optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 Lb/2],[F*2 1.42],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));

x
fval
F = x(1);
d = x(2);
Lb = x(3);


% Lb_min=0.5;%1.2;%
% Lb_max=1.4;%2.2;%2.8;
% nptLb = 19;%31;%41;
% 
% disp('start')
% tps0=now;
% datestr(tps0)
% tic
% for i=1:nptLb
%     i
%     Lb = Lb_min+(i-1)*(Lb_max-Lb_min)/(nptLb-1);
%     %[x,fval,exitflag] = fminsearchbnd(@null,x0,[0 0],[1-2*d*tan(pente)/Lb (1-F)*Lb/(2*tan(pente))],optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%     [x,fval,exitflag] = fminsearchbnd(@null,x0,[0 0],[2*F 2*d],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%     F = x(1);
%     d = x(2);
%     x0 = x
%     xx(i) = Lb;
%     yy1(i) = F;
%     yy2(i) = d;
%     yy3(i) = mean(nul_res_sp_b+Nghost);
%     yy4(i) = mean(nul_res_sp_b+NARGghost);
%     yy5(i) = mean(nul_res_sp_b);
%     yy6(i) = mean(tmp3);
%     yy7(i) = mean(Tin);
%     yy8(i) = mean(TARGtot);
%     %toc
%     tpsfin = tps0 + nptLb/i*(now-tps0);
%     disp('Fin prévue (date+heure) :')
%     datestr(tpsfin)
% end
% disp('end')
% datestr(now)
% toc
% 
% 
% 
% % Figures
% % -------
% 
% close all
% newFig
% set(gca,'XLim',[Lb_min Lb_max])%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% xlabel('Period  $\Lambda$  $(\mu m)$','FontSize',fsz)
% %ylabel('Fill factor','FontSize',fsz)
% [ax,h1,h2] = plotyy(xx,yy1,xx,yy2);
% set(ax(2),'xtick',[])
% set(ax(1),'XLim',[Lb_min Lb_max],'ycolor','k')
% set(ax(2),'XLim',[Lb_min Lb_max],'ycolor','k')
% set(h1,'color',myblue,'linewidth',lwz)
% set(h2,'color',mygreen,'linewidth',lwz)
% ylabel(ax(1),'Fill factor','FontSize',fsz*1.2,'color','k')
% ylabel(ax(2),'Grating depth  $(\mu m)$','FontSize',fsz*1.2,'color','k')
% leg=legend(' Fill',' Depth');
% set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
% tick2latex
% print('-depsc2',sprintf('Lb_F+d.eps'), '-r300');
% 
% 
% 
% newFig
% set(gca,'XLim',[Lb_min Lb_max])%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% xlabel('Period  $\Lambda$  $(\mu m)$','FontSize',fsz)
% ylabel('Null depth','FontSize',fsz)
% plot(xx,yy3,'-','color',myred,'linewidth',lwz)
% plot(xx,yy4,'-','color',mygreen,'linewidth',2*lwz)
% plot(xx,yy5,'--','color',myblue,'linewidth',1.2*lwz)
% leg=legend(' without ARG',' with ARG',' no ghost');
% set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
% tick2latex
% print('-depsc2',sprintf('Lb_Null.eps'), '-r300');
% 
% newFig
% set(gca,'XLim',[Lb_min Lb_max])%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% xlabel('Period  $\Lambda$  $(\mu m)$','FontSize',fsz)
% ylabel('Phase shift $\Phi_{\rm{TE-TM}}$','FontSize',fsz)
% plot(xx,xx.*0+pi,'--','color',myblue,'linewidth',lwz)
% plot(xx,yy6,'-','color',mygreen,'linewidth',2*lwz)
% tick2latex
% print('-depsc2',sprintf('Lb_Phase.eps'), '-r300');
% 
% newFig
% set(gca,'XLim',[Lb_min Lb_max])%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% set(gca,'ylim',[.6 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% xlabel('Period  $\Lambda$  $(\mu m)$','FontSize',fsz)
% ylabel('Transmittance','FontSize',fsz)
% plot(xx,yy7,'--','color',myred,'linewidth',lwz)
% plot(xx,yy8,'-','color',mygreen,'linewidth',2*lwz)%tit=title('L-band AGPM with ARG');
% %set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
% leg=legend(' $T_{\rm{out}} = T_{\rm{ARG}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' Total trans.');
% set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','south')
% leg=legend(' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' Total trans.');
% set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
% tick2latex
% print('-depsc2',sprintf('Lb_Trans.eps'), '-r300');
% 
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
y1=nul_res_sp_b+Nghost;
y2=nul_res_sp_b+NARGghost;
y3=nul_res_sp_b;
plot(lb_t,y1,'-','color',myred,'linewidth',lwz)
plot(lb_t,y2,'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,y3,'--','color',myblue,'linewidth',1.2*lwz)
leg=legend(' without ARG',' with ARG',' no ghost');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_Null.eps'), '-r300');
% mean(y1(lb_t<=4.1 | lb_t>=4.6))
% mean(y2(lb_t<=4.1 | lb_t>=4.6))
% mean(y3(lb_t<=4.1 | lb_t>=4.6))
% mean(y1(lb_t<=4.1))
% mean(y2(lb_t<=4.1))
% mean(y3(lb_t<=4.1))
% mean(y1(lb_t>=4.6))
% mean(y2(lb_t>=4.6))
% mean(y3(lb_t>=4.6))



% Phase TE-TM
% -----------
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Phase shift $\Phi_{\rm{TE-TM}}$','FontSize',fsz)
%plot(lb_t,nul_res_sp_b+Nghost,'-','color',myred,'linewidth',lwz)
plot(lb_t,lb_t.*0+pi,'--','color',myblue,'linewidth',lwz)
plot(lb_t,tmp3,'-','color',mygreen,'linewidth',2*lwz)
%str = sprintf('Outer L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Full L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Central L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%title(str,'FontSize',fsz)
%leg=legend(' Without ARG',' With ARG',' Optimal');
%set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_PhaseTE-TM.eps'), '-r300');



% Total Trans WITH ARG
% --------------------
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'ylim',[.7 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
% plot(lb_t,pchip(lb_t,TARG,lb_t),'--','color',myblue,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,Tin,lb_t),'--','color',myred,'linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp2,lb_t),'k+:','linewidth',lwz)
% plot(lb_t,pchip(lb_t,tmp1,lb_t),'k.:','linewidth',lwz)
% plot(lb_t,pchip(lb_t,TARGtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,TARG,'-','color',myblue,'linewidth',lwz)
plot(lb_t,Tin,'-','color',myred,'linewidth',lwz)
plot(lb_t,tmp2,'k-.','linewidth',lwz)
plot(lb_t,tmp1,'k--','linewidth',lwz)
plot(lb_t,TARGtot,'-','color',mygreen,'linewidth',2*lwz)%tit=title('L-band AGPM with ARG');
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{out}} = T_{\rm{ARG}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_TtotARG.eps'), '-r300');




break


% Total Trans WITHOUT ARG
% -----------------------
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
tit=title('L-band AGPM w/o ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{out}} = \frac{4n}{(n+1)^2}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_Ttot.eps'), '-r300');



% Total Refl WITH ARG
% -------------------
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
set(gca,'ylim',[0 .2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
plot(lb_t,pchip(lb_t,RARGtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,RARG,lb_t),'--','color',myblue,'linewidth',lwz)
tit=title('L-band AGPM with ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' Total refl.',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$',' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',' $R_{\rm{out}} = R_{\rm{ARG}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_RtotARG.eps'), '-r300');



% Total Refl WITHOUT ARG
% ----------------------
figure('name','relf')
set(gcf,'color',[1 1 1])
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'ylim',[0 .4])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Reflectance','FontSize',fsz)
plot(lb_t,pchip(lb_t,Rtot,lb_t),'-','color',mygreen,'linewidth',2*lwz)
plot(lb_t,pchip(lb_t,R0,lb_t),'--','color',myblue,'linewidth',lwz)
plot(lb_t,pchip(lb_t,Rin,lb_t),'--','color',myred,'linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp4,lb_t),'k.:','linewidth',lwz)
plot(lb_t,pchip(lb_t,tmp5,lb_t),'k+:','linewidth',lwz)
tit=title('L-band AGPM w/o ARG');
set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' Total refl.',' $R_{\rm{out}} = \frac{(n-1)^2}{(n+1)^2}$',' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
tick2latex
print('-depsc2',sprintf('agpm_L_Rtot.eps'), '-r300');





NullNO=mean(nul_res_sp_b)
NullYES=mean(nul_res_sp_b+Nghost);
NullYESARG=mean(nul_res_sp_b+NARGghost)


Nghostmean=mean(Nghost)
NARGghostmean=mean(NARGghost);







