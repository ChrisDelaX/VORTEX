%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Lib�ration 
%----------
save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost bandAR nul_ghost nul_ghost_ARG nul_ideal
global tmp1p1 tmp1m1 tmp2p1 tmp2m1


%Allocation m�moire
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
%13: profil trap�zoidal simple
%   -> L couches
%14: profil trap�zoidal avec couche antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2 couches
%15: profil trap�zoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trap�zoidal avec couche antireflet continue
%   -> L+1 couches
%17: profil trap�zoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trap�zoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2+1 couches
%19: profil trap�zoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trap�zoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trap�zoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trap�zoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2+1 couches
%23: profil trap�zoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trap�zoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)
p=13;%1;%
%Nombre de couches de discr�tisation du profil non rectangulaire en sus des milieux ext�rieurs
L=16;%32;%25;%50;%


%Choix mat�riau --> Permittivit�s
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
%Permittivit�s milieu ext�rieur incident
EI_choice=1;%3;%
%Permittivit�s milieu ext�rieur �mergent
EIII_choice=3;%1;%
%Permittivit�s du r�seau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;
%Permittivit�s de la couche antireflet
E_AR_choice=1;%14;
E_AR1_choice=1;
E_AR2_choice=1;
%Permittivit�s de la couche d'arret
E_arret_choice=1;


%Onde incidente
%--------------
%Angle incidence non conique
theta=deg2rad(0); theta_min=deg2rad(30); theta_max=deg2rad(50);
%Angle incidence conique
phi=deg2rad(0); phi_min=deg2rad(0); phi_max=deg2rad(0);
%Angle polarisation
psi=deg2rad(45); psi_min=deg2rad(0); psi_max=deg2rad(0);


% Absorption
% ----------
nlb=11;%37;%81;%116;%
lb_min=2;%1.9;%2;%9;%11;%3.5;%3.5;%4.6;%8;%
lb_max=13.5;%2.4;%2.5;%13.2;%4.1;%5.0;%11.3;%
lb_t=lb_min:(lb_max-lb_min)/(nlb-1):lb_max;
bandx = [lb_min lb_max];
K0 = 2;K1 = 4;
N0 = 8;N1 = 13.5;
%lb_t = [K0:(K1-K0)/(nlb-1):K1 N0:(N1-N0)/(nlb-1)/2:N1];
lb_t = [K0:(K1-K0)/(nlb-1):K1 N0:(N1-N0)/(nlb-1):N1];
TH_diam_new

% Anti-reflection grating
% -----------------------
lar_min = 10;
lar_max = 12.5;
bandAR = [lar_min (lar_max+lar_min)/2 lar_max];


%Param�tres du r�seau
%--------------------
% p�riode --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
Lb=3.6474;%3.5891;%4.6;%1.4;%2.5;%1.42;%1.89;%3.32;%4.58;%
F=.5397;%.5362;%.4826;%.4;%1.0/Lb;%0.37;%0.4506;%0.4142;%0.4930;%0.4449;%
d=17.5617;%16.9154;%15.83;%13.8;%13.5;%8.4;%4.76;%5.2226;%6.0735;%15.7734;%16.5493;%
%x0=[F d];
x0=[Lb F d];
%x0=[];
pente_deg=2.45;%2.75;%5.0;%3.2;%2.5;%3.25;%2.75;%
pente=deg2rad(pente_deg);

%COND LIMIT
Fmax = 1-2*d*tan(pente)/Lb
dmax = (1-F)*Lb/(2*tan(pente))
%x0=[Fmax dmax];

%
% Calcul rcwa
% -----------
%[fval1,fval2] = null();
%[x,fval,exitflag] = fminsearchbnd(@null,x0,Lb/2,Lb,optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,F/2,F*2,optimset('MaxIter',25,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,d/2,d*2,optimset('MaxIter',25,'Display','iter','TolX',1e-5,'TolFun',1e-6));
%[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));
[x,fval,exitflag] = fminsearchbnd(@null,x0,[Lb/2 F/2 d/2],[Lb*1.2 F*2 d*2],optimset('MaxIter',15,'Display','iter','TolX',1e-5,'TolFun',1e-6));

%x
%fval
%save('data.txt','lb_t','tmp3','Ttot','TARGtot','nul_ghost','nul_ghost_ARG','nul_ideal','-ascii')
return

% 2D MAP
% ------
F_min = .48;F_max = .56;
d_min = 14;d_max = 17.5;
npts_F = 25;
npts_d = 25;

tic
depart_tmp=now;
depart=datestr(depart_tmp)
for j=1:npts_F
    j
    F = F_min+(F_max-F_min)*(j-1)/(npts_F-1);
    Xparam(j)=F;
    for k=1:npts_d
        d = d_min+(d_max-d_min)*(k-1)/(npts_d-1);
        Yparam(k)=d;
        [fval1,fval2] = null();
        if imag(fval1) == 0
            nuldpt(k,j) = mean(fval1);
        else
            nuldpt(k,j) = 0;
        end
        if j==1 && k==1
            fin=datestr(depart_tmp+toc/3600/24*npts_F*npts_d)
        end
    end
end
fmulti = nuldpt;    
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
xlabel('Fill factor','FontSize',fsz)
ylabel('Grating depth $(\mu m)$','FontSize',fsz)

Zparam = (fmulti');
colorbarLim = [0 1];
Zparam(Zparam>colorbarLim(2)) = colorbarLim(2);
hS = surfc(Xparam,Yparam,Zparam);
%shading('interp') %flat, faceted
set(hS,'EdgeColor', 'none')
hC = colorbar;
set(hC,'box','on','linewidth',lwz/2)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(hC,'YLim',colorbarLim)
print('-depsc2',sprintf('map2D.eps'), '-r300');

return

%contours
[C01,C02] = contour(Xparam,Yparam,Zparam,[-3.8:.2:-3.0],'k-');
[C03,C04] = contour(Xparam,Yparam,Zparam,[-2.7:.5:-1.7],'k--');
%[C05,C06] = contour(Xparam,Yparam,Zparam,[-1.6:.4:-0],'k-.');
clabel(C01,C02,'FontSize',fsz);
clabel(C03,C04,'FontSize',fsz);
clabel(C05,C06,'FontSize',fsz);


return

% Figures
% -------


% Null Depth + GHOST
% ------------------
close all
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[1e-4 1e-1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Peak raw attenuation (null depth)','FontSize',fsz)
lb_N = lb_t(nlb+1:end)
set(gca,'XLim',[lb_N(1) lb_N(end)])
plot(lb_N,nul_res_sp_b(nlb+1:end),'-','color',myblue,'linewidth',1.2*lwz)
%plot(lb_t,nul_ghost,'-','color',myred,'linewidth',lwz)
%plot(lb_t,nul_ghost_ARG,'-','color',mygreen,'linewidth',2*lwz)
%plot(lb_t,nul_ideal,'-','color',myblue,'linewidth',1.2*lwz)
%plot(lb_t,polyval(polyfit(lb_t,nul_ghost_ARG,2),lb_t),'-','color',mygreen,'linewidth',2*lwz)
%str = sprintf('Outer L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('\Lambda = %3.3g \mu m, h = %3.3g \mu m,  F = %.2g',Lb,d,F);
%str = sprintf('Central L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%title(str,'FontSize',fsz)
%leg=legend(' without ARG',' with ARG',' no ghost');
%leg=legend(' w/o ARG',' w/ ARG');
%set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','northeast')
%tick2latex
%save
%print('-depsc2',sprintf('agpm_L_outer.eps'), '-r300');
print('-depsc2',sprintf('agpm_Nband_Null.eps'), '-r300');
%print('-dpng',sprintf('agpm_L_brun.png'), '-r300');
%break


return

% Orders
% ------

newFig
%set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[0 2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Diffraction orders','FontSize',fsz)
Tin = (tmp1(1:nlb)+tmp2(1:nlb))./2
Tinp1 = (tmp1p1(1:nlb)+tmp2p1(1:nlb))./2
Tinm1 = (tmp1m1(1:nlb)+tmp2m1(1:nlb))./2
lb_K = lb_t(1:nlb)
set(gca,'XLim',[lb_K(1) lb_K(end)])
plot(lb_K,Tin,':k','linewidth',lwz)
plot(lb_K,Tinp1,'-','color',mygreen,'linewidth',2*lwz)
plot(lb_K,Tinm1,'-','color',myred,'linewidth',2*lwz)
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $m = 0$',' $m = +1$',' $m = -1$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
%tick2latex
print('-depsc2',sprintf('agpm_Nband_orders.eps'), '-r300');

return

tmp3 = (mod((tmp3+pi),2*pi)-pi)./pi;

for i=2:nlb-1
    if abs(tmp1(i)-tmp1(i-1))>.04
        tmp1(i)=(tmp1(i-1)+tmp1(i+1))./2;
    end
    if abs(tmp2(i)-tmp2(i-1))>.04
        tmp2(i)=(tmp2(i-1)+tmp2(i+1))./2;
    end
    if abs(tmp3(i)-tmp3(i-1))>.04
        tmp3(i)=(tmp3(i-1)+tmp3(i+1))./2;
    end
end


% Phase TE-TM
% -----------
newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
set(gca,'ylim',[-.8 .8])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
%set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Phase shift $\Phi_{\rm{TE-TM}} \: \: \: (\pi \rm{rad})$','FontSize',fsz)
%plot(lb_t,lb_t.*0+1,'--','color',myblue,'linewidth',lwz)
%plot(lb_t,tmp3./pi,'-','color',mygreen,'linewidth',2*lwz)
lb_K = lb_t(1:nlb)
set(gca,'XLim',[lb_K(1) lb_K(end)])
phi = (mod(tmp3(1:nlb)+pi,2*pi)-pi)./pi
plot(lb_K,phi,'-','color',mygreen,'linewidth',2*lwz)
%plot(lb_t,polyval(polyfit(lb_t,tmp3,2),lb_t),'-','color',mygreen,'linewidth',2*lwz)
%str = sprintf('Outer L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Full L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Central L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%title(str,'FontSize',fsz)
%leg=legend(' Without ARG',' With ARG',' Optimal');
%set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
%tick2latex
%save
%print('-depsc2',sprintf('agpm_L_outer.eps'), '-r300');
print('-depsc2',sprintf('agpm_Nband_PhaseTE-TM.eps'), '-r300');
%print('-depsc2',sprintf('agpm_L_central.eps'), '-r300');
% break


% Total Trans WITH ARG
% --------------------

newFig
set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
%set(gca,'ylim',[0 2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
xlabel('Wavelength  $\lambda$  $(\mu m)$','FontSize',fsz)
ylabel('Transmittance','FontSize',fsz)
plot(lb_t,tmp2,'k+:','linewidth',lwz)
plot(lb_t,tmp1,'k.:','linewidth',lwz)
plot(lb_t,Tin,'--','color',myred,'linewidth',lwz)
plot(lb_t,TARG,'--','color',myblue,'linewidth',lwz)
plot(lb_t,TARGtot,'-','color',mygreen,'linewidth',2*lwz)
%tm = polyval(polyfit(lb_t,tmp2,2),lb_t);
%te = polyval(polyfit(lb_t,tmp1,2),lb_t);
%plot(lb_t,tm,'k-.','linewidth',lwz)
%plot(lb_t,te,'k--','linewidth',lwz)
%plot(lb_t,tm./te,'-','color',mygreen,'linewidth',2*lwz)%tit=title('L-band AGPM with ARG');
%set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',' $T_{\rm{out}} = T_{\rm{ARG}}$',' Total trans.');
%leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{TM}}/T_{\rm{TE}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','best')
%tick2latex
print('-depsc2',sprintf('agpm_Nband_TtotARG.eps'), '-r300');

return



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
%tick2latex
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
%tick2latex
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
%tick2latex
print('-depsc2',sprintf('agpm_L_Rtot.eps'), '-r300');





NullNO=mean(nul_res_sp_b)
NullYES=mean(nul_res_sp_b+Nghost);
NullYESARG=mean(nul_res_sp_b+NARGghost)


Nghostmean=mean(Nghost)
NARGghostmean=mean(NARGghost);







