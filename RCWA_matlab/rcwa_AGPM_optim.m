%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Libï¿½ration 
%----------
save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 
global d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 
global E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice 
global E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente 
global fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost 
global bandAR bandTAR nul_ghost nul_ghost_ARG nul_ideal
global tmp1p1 tmp1m1 tmp2p1 tmp2m1 optim


%Allocation mï¿½moire
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
%13: profil trapï¿½zoidal simple
%   -> L couches
%14: profil trapï¿½zoidal avec couche antireflet dans les creux, sur les bosses et ...
%   adhï¿½rent aux parois obliques
%   -> L+2 couches
%15: profil trapï¿½zoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trapï¿½zoidal avec couche antireflet continue
%   -> L+1 couches
%17: profil trapï¿½zoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trapï¿½zoidal avec couche d'arret continue et antireflet dans les creux, 
%   sur les bosses et adhï¿½rent aux parois obliques
%   -> L+2+1 couches
%19: profil trapï¿½zoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trapï¿½zoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trapï¿½zoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trapï¿½zoidal avec couche d'arret discontinue et antireflet dans les creux, 
%   sur les bosses et adhï¿½rent aux parois obliques
%   -> L+2+1 couches
%23: profil trapï¿½zoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trapï¿½zoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)
p=13;%1;%
%Nombre de couches de discrï¿½tisation du profil non rectangulaire en sus des milieux extï¿½rieurs
L=25;%16;%32;%50;%


%Choix matï¿½riau --> Permittivitï¿½s
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
%Permittivitï¿½s milieu extï¿½rieur incident
EI_choice=1;%3;%
%Permittivitï¿½s milieu extï¿½rieur ï¿½mergent
EIII_choice=3;%1;%
%Permittivitï¿½s du rï¿½seau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;
%Permittivitï¿½s de la couche antireflet
E_AR_choice=1;%14;
E_AR1_choice=1;
E_AR2_choice=1;
%Permittivitï¿½s de la couche d'arret
E_arret_choice=1;


%Onde incidente
%--------------
%Angle incidence non conique
theta=deg2rad(0); theta_min=deg2rad(30); theta_max=deg2rad(50);
%Angle incidence conique
phi=deg2rad(0); phi_min=deg2rad(0); phi_max=deg2rad(0);
%Angle polarisation
psi=deg2rad(45); psi_min=deg2rad(0); psi_max=deg2rad(0);


% Select parameters
% ------------------
name = 'AGPM-METIS';
bands = string({'L','M','N1','N2'});%
nlb = 81;
optim = 2; % 0:d; 1:F,d; 2:nothing
ARtrans = [.993 .999 .995];

for band = bands

switch band
    case 'L'
        lam = create_param('Wavelength', 'µm', 3.5, 1.2/3.5, nlb);
        ARlam = create_param('AR wavelength', 'µm', lam.val, 0.2, 3);
        per = create_param('Period', 'µm', 1.21, 0, 1);
        lw = create_param('Line width', 'µm', 0.65, 0, 1);
        dep = create_param('Depth', 'µm', 5.53, 0, 1);
        sw = create_param('Sidewall angle', 'deg', 2.45, 0, 1);
    case 'M'
        lam = create_param('Wavelength', 'µm', 4.6, 1.4/4.6, nlb);
        ARlam = create_param('AR wavelength', 'µm', lam.val, 0.2, 3);
        per = create_param('Period', 'µm', 1.63, 0, 1);
        lw = create_param('Line width', 'µm', 0.83, 0, 1);
        dep = create_param('Depth', 'µm', 6.55, 0, 1);
        sw = create_param('Sidewall angle', 'deg', 2.45, 0, 1);
    case 'N1'
        lam = create_param('Wavelength', 'µm', 9.25, 2.5/9.25, nlb);
        ARlam = create_param('AR wavelength', 'µm', lam.val, 0.2, 3);
        per = create_param('Period', 'µm', 3.36, 0, 1);
        lw = create_param('Line width', 'µm', 1.64, 0, 1);
        dep = create_param('Depth', 'µm', 12.57, 0, 1);
        sw = create_param('Sidewall angle', 'deg', 2.45, 0, 1);
    case 'N2'
        lam = create_param('Wavelength', 'µm', 11.75, 3.5/11.75, nlb);
        ARlam = create_param('AR wavelength', 'µm', lam.val, 0.2, 3);
        per = create_param('Period', 'µm', 4.2, 0, 1);
        lw = create_param('Line width', 'µm', 2.11, 0, 1);
        dep = create_param('Depth', 'µm', 16.57, 0, 1);
        sw = create_param('Sidewall angle', 'deg', 2.45, 0, 1);
end

% populate global variables
lb_t = lam.vals;
bandAR = ARlam.vals;
bandTAR = ARtrans;
Lb = per.val;
F = lw.val/Lb;
d = dep.val;
pente_deg = sw.val;
pente = deg2rad(pente_deg);    
        
% material absorption
TH_diam_new

%% Start optimization
% -------------------
ncurves = 1;
lw.step = 0.01;
for ii=1:ncurves
lw.temp = lw.val + (ii-(ncurves+1)/2)*lw.step;
F = lw.temp/Lb;
d = dep.val;
ii
switch optim
    case 0
        x0 = [d]
    case 1
        x0 = [F, d]
    case 2
        x0 = [];
end

%COND LIMIT
Fmax = 1-2*d*tan(pente)/Lb;
dmax = (1-F)*Lb/(2*tan(pente));
%x0=[Fmax dmax];
if ii == 1
    Fmax
    dmax
end

%
% Calcul rcwa
% -----------
switch optim
    case 0
        [x,fval,exitflag] = fminsearchbnd(@null,x0,d/2,d*2,optimset('MaxIter',15,...
            'Display','iter','TolX',1e-5,'TolFun',1e-6));
        d = x(1)
    case 1
        [x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',15,...
            'Display','iter','TolX',1e-5,'TolFun',1e-6));
        F = x(1)
        d = x(2)
    case 2
        [fval,fval2] = null();
end
fval

%save('data.txt','lb_t','tmp3','Ttot','TARGtot','nul_ghost','nul_ghost_ARG','nul_ideal','-ascii')
nulls(ii,:)=nul_ideal;
lws(ii)=F*Lb;
deps(ii)=d;
if ii <= 3
    opt_ghost = nul_ghost;
    opt_ARG = nul_ghost_ARG;
    opt_null = nul_ideal;
    opt_lw = F*Lb;
    opt_d = d;
end
end

%% Null Depth + GHOST
% ------------------
close all
newFig
xlabel(sprintf('%s (%s)',lam.name,lam.units),'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz*1.2)
ylabel('Peak raw contrast (null depth)','Fontname',fnz,'FontWeight',fwz,'FontSize',fsz*1.2)
title(sprintf('%s band [%3.2f-%3.2f %s], %s = %3.2f %s,\n %s = %3.2f %s, %s = %3.2f %s, %s = %3.2f %s',...
    band, lam.range(1), lam.range(2), lam.units, per.name, per.val, per.units,...
    lw.name, opt_lw, lw.units, dep.name, opt_d, dep.units, sw.name, sw.val, sw.units), ...
    'Fontname',fnz,'FontSize',fsz,'FontWeight','bold');%,'horizontalalignment','left')
set(gca,'XLim',lam.range)%,'xtick',bandtick)
set(gca,'YLim',[1e-4 0.01])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
%plot(lb_t,opt_ghost,':','color',myred,'linewidth',lwz)
plot(lb_t,opt_ARG,'-.','color',myblue,'linewidth',lwz*1.2)
plot(lb_t,opt_null,'-','color',mygreen,'linewidth',lwz*1.2)
leg=legend('  ghost included','  ideal (no ghost)');
%leg=legend('  without ARG','  with ARG','  ideal (no ghost)');
set(leg,'box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz*1.2,'location','northeast')
%tick2latex
print('-dpng',sprintf('%s-%s_1D.png',name,band), '-r300');
fitswrite([lb_t; opt_ARG; opt_null],sprintf('%s-%s_1D.fits',name,band))
end
return

%% Multi linewidths
% -----------------

close all
newFig
xlabel(sprintf('%s (%s)',lam.name,lam.units),'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz*1.2)
ylabel('Peak raw contrast (null depth)','Fontname',fnz,'FontWeight',fwz,'FontSize',fsz*1.2)
title(sprintf('%s band [%3.2f-%3.2f %s], %s = %3.2f %s, %s = %3.2f %s', band, ...
    lam.range(1), lam.range(2), lam.units, per.name, per.val, per.units,...
    sw.name, sw.val, sw.units), 'Fontname',fnz,'FontSize',fsz,'FontWeight','bold');%,'horizontalalignment','left')
set(gca,'XLim',lam.range)%,'xtick',bandtick)
set(gca,'YLim',[1e-4 0.02])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
plot(lb_t,nulls(1,:),':','color','m','linewidth',lwz)
plot(lb_t,nulls(2,:),'-.','color',myred,'linewidth',lwz)
plot(lb_t,nulls(3,:),'-','color',mygreen,'linewidth',lwz*1.2)
plot(lb_t,nulls(4,:),'-.','color',myblue,'linewidth',lwz)
plot(lb_t,nulls(5,:),':','color','c','linewidth',lwz)
str1 = sprintf('  %s = %3.2f %s, %s = %3.2f %s',lw.name,lws(1),lw.units,dep.name,deps(1),dep.units);
str2 = sprintf('  %s = %3.2f %s, %s = %3.2f %s',lw.name,lws(2),lw.units,dep.name,deps(2),dep.units);
str3 = sprintf('  %s = %3.2f %s, %s = %3.2f %s',lw.name,lws(3),lw.units,dep.name,deps(3),dep.units);
str4 = sprintf('  %s = %3.2f %s, %s = %3.2f %s',lw.name,lws(4),lw.units,dep.name,deps(4),dep.units);
str5 = sprintf('  %s = %3.2f %s, %s = %3.2f %s',lw.name,lws(5),lw.units,dep.name,deps(5),dep.units);
leg=legend(str1,str2,str3,str4,str5);
set(leg,'box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',...
    fwz,'FontSize',fsz*1.2,'location','northeast')
%tick2latex
print('-dpng',sprintf('%s-%s_1D_multi.png',name,band), '-r300');




%% Orders
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
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz,'location','best')
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
%str = sprintf('Outer L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: ...
%   h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Full L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: ...
%   h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%str = sprintf('Central L band $$\\rightarrow \\Lambda = %3.3g \\mu m, \\:\\: ...
%   h = %3.3g \\mu m \\rightarrow F = %.2g $$',Lb,d,F);
%title(str,'FontSize',fsz)
%leg=legend(' Without ARG',' With ARG',' Optimal');
%set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
%   'FontSize',fsz,'location','best')
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
leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',...
    ' $T_{\rm{out}} = T_{\rm{ARG}}$',' Total trans.');
%leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{TM}}/T_{\rm{TE}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz,'location','best')
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
leg=legend(' $T_{\rm{TM}}$',' $T_{\rm{TE}}$',' $T_{\rm{in}} = \frac{T_{\rm{TE\;}}+T_{\rm{TM}}}{2}$',...
    ' $T_{\rm{out}} = \frac{4n}{(n+1)^2}$',' Total trans.');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz,'location','best')
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
leg=legend(' Total refl.',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$',...
    ' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',...
    ' $R_{\rm{out}} = R_{\rm{ARG}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz,'location','best')
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
leg=legend(' Total refl.',' $R_{\rm{out}} = \frac{(n-1)^2}{(n+1)^2}$',...
    ' $R_{\rm{in}} = \frac{R_{\rm{TE\;}}+R_{\rm{TM}}}{2}$',' $R_{\rm{TE}}$',' $R_{\rm{TM}}$');
set(leg,'interpreter','latex','box','on','linewidth',lwz,'Fontname',fnz,'FontWeight',fwz,...
    'FontSize',fsz,'location','best')
%tick2latex
print('-depsc2',sprintf('agpm_L_Rtot.eps'), '-r300');





NullNO=mean(nul_res_sp_b)
NullYES=mean(nul_res_sp_b+Nghost);
NullYESARG=mean(nul_res_sp_b+NARGghost)


Nghostmean=mean(Nghost)
NARGghostmean=mean(NARGghost);







