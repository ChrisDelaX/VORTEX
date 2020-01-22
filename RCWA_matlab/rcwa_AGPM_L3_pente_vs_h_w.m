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
warning off MATLAB:polyfit
% figures parameters
fnz = 'Arial'; % fontname
fsz = 26; % fontsize
fwz = 'normal';%'Bold'; % fontweight
lwz = 2;  % line width


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

p=31;%13;%1;%3;%

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=25;%50;%8;%
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
Lb=1.42;%4.6;%
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
lb_min=3.25;%3.5;%
lb_max=4.25;%4.1;%
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


%Optimisation
%------------


% F=0.4;
% d=4.3;%5.05;
% x0=[d];
% [x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
% x

% valeurs correctes: 5.1979    0.4651
% valeurs incorrectes: 4.3031    0.4026

disp('start')
datestr(now)
tic

x0=[];

% interpolation linéaire de la pente alpha
x1= 4.7;
y1= 2.9;
x2= 5.57;
y2= 0.85*rad2deg(atan((1-0.41)*1.42/2/x2));


pen_a = (y2-y1)/(x2-x1);
pen_b = y1 - pen_a*x1;

% Filters
% -------

load filters.mat
x1 = filter_wide(:,1);
y1 = filter_wide(:,2);
p_wide = fit(x1,y1,'smoothingspline');


% AGPM-L3
% -------

L=25;%100;%2;%
F=0.4576;%0.4404;%
d=6;%5.8;%
%F=F+d/Lb*tan(rad(3.25));
pen_max = rad2deg(atan((1-F)*Lb/2/(d*5/6)));


% figures optim
% -------------

pente=rad(3.95);
d=6;
F=0.65/1.42;

d_min = 5.2;%3.25;%
d_max = 6.2;
w_min = 0.62;
w_max = 0.74;
npt_d=10;%50;%
npt_w=10;%50;%
for i=1:npt_d
    i
    xd(i) = d_min + (d_max-d_min)*(i-1)/(npt_d-1);
    d=xd(i);
    pente=rad(3.9);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y1d(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
    pente=rad(3.8);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y2d(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
    pente=rad(4.0);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y3d(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
end
for i=1:npt_w
    i
    xw(i) = w_min + (w_max-w_min)*(i-1)/(npt_w-1);
    F=xw(i)/Lb;
    pente=rad(3.95);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y1w(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
    pente=rad(3.80);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y2w(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
    pente=rad(4.10);
    [fval1,fval2] = null(x0);
    nul_tmp3 = nul_res_sp_b;
    y3w(i) = sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
end
    
close all
% h +- delta h
figure
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'YScale','log')
set(gca,'XLim',[d_min d_max],'ylim',[.00001 .1])
%set(gca,'xtick',[1 200 400 600 800 1000])
%title('Schematic profile','FontSize',fsz*1.2)
xlabel('Total height h (µm)','FontSize',fsz*1.1)
ylabel('Raw null depth','FontSize',fsz*1.1)
plot(xd,y1d,'color',[0.2 0.6 0.2],'linewidth',lwz*2)
plot(xd,y2d,'b--','linewidth',lwz*2)
plot(xd,y3d,'r-.','linewidth',lwz*2)
l=legend(' \alpha=3.9°',' \alpha-\Delta\alpha',' \alpha+\Delta\alpha');
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)

print('-depsc2',sprintf('deltah.eps'), '-r300');

%break

% w +- delta w
figure
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'YScale','log')
set(gca,'XLim',[w_min w_max],'ylim',[.00001 .1])
%set(gca,'xtick',[1 200 400 600 800 1000])
%title('Schematic profile','FontSize',fsz*1.2)
xlabel('Upper width w (µm)','FontSize',fsz*1.1)
ylabel('Raw null depth','FontSize',fsz*1.1)
plot(xw,y1w,'color',[0.2 0.6 0.2],'linewidth',lwz*2)
plot(xw,y2w,'b--','linewidth',lwz*2)
plot(xw,y3w,'r-.','linewidth',lwz*2)
l=legend(' \alpha=3.9°',' \alpha-\Delta\alpha',' \alpha+\Delta\alpha');
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)

print('-depsc2',sprintf('deltaw.eps'), '-r300');


break


% calcul réjection, cas réalist
pen_pontus = rad2deg(atan(((1-F)*Lb-0.080)/2/(d*5/6)))
penL3=pen_pontus;%4.1;%3.25;%0;%2.9;%pen_max;%
pente=rad(penL3);
[fval1,fval2] = null(x0);
nul_tmp3 = nul_res_sp_b;
FL3=F;
hL3=d;
%RL3 = 1/mean(nul_tmp3)%*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
NL3=sum(nul_tmp3.*p_wide(lb_t)')/sum(p_wide(lb_t));
RL3_filter = 1/NL3

break

% cas optimist
pen_opt = rad2deg(atan(((1-F)*Lb-0.060)/2/(d*5/6)))
pente=rad(pen_opt);
[fval1,fval2] = null(x0);
nul_tmp3_opt = nul_res_sp_b;
NL3_opt=sum(nul_tmp3_opt.*p_wide(lb_t)')/sum(p_wide(lb_t));
RL3_filter_opt = 1/NL3_opt


% cas pesimist
pen_pes = rad2deg(atan(((1-F)*Lb-0.100)/2/(d*5/6)))
pente=rad(pen_pes);
[fval1,fval2] = null(x0);
nul_tmp3_pes = nul_res_sp_b;
NL3_pes=sum(nul_tmp3_pes.*p_wide(lb_t)')/sum(p_wide(lb_t));
RL3_filter_pes = 1/NL3_pes


disp('end')
datestr(now)
toc



% Figures
% *******
% *******


% Profile
% *******

%close all
figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
%set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
prof_fig=flipdim(prof,1);

subplot(1,2,1)
imagesc(prof_fig)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',[1 1024],'ylim',[1 L])
set(gca,'xtick',[1 200 400 600 800 1000])
title('Schematic profile','FontSize',fsz*1.2)
xlabel('Sampled period (1024 samples)','FontSize',fsz*1.1)
ylabel(sprintf('Layer number (L=%3.0f)',L),'FontSize',fsz*1.1)
if L==2
    set(gca,'ytick',[1])
    ylabel(sprintf('One layer (L=%3.0f)',L-1),'FontSize',fsz*1.1)
end


subplot(1,2,2)
imagesc(prof_fig)
set(gca,'PlotBoxAspectRatio',[1 6/1.42 1])
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',[1 1024],'ylim',[1 L])
set(gca,'xtick',[1 1024])
title('Scaled profile','FontSize',fsz*1.2)
xlabel('Sampled period','FontSize',fsz*1.1)
ylabel(sprintf('Layer number (L=%3.0f)',L),'FontSize',fsz*1.1)
if L==2
    set(gca,'ytick',[1])
    ylabel(sprintf('One layer (L=%3.0f)',L-1),'FontSize',fsz*1.1)
end

%colormap(Summer)
colormap(Copper)
%text(400,1.6,,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
gtext({sprintf('slope=%3.2f°',penL3),sprintf('top width=%3.2fµm',F*Lb),sprintf('depth=%3.2fµm',d),' ',sprintf('rejection=%3.2f',RL3_filter),sprintf('null depth=%3.4f',NL3)},'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)

% print('-dpng',sprintf('prof_rect_equiv.png'), '-r300');
% print('-depsc2',sprintf('prof_rect_equiv.eps'), '-r300');

% print('-dpng',sprintf('prof_trap.png'), '-r300');
% print('-depsc2',sprintf('prof_trap.eps'), '-r300');

% print('-dpng',sprintf('prof_trap+tri100.png'), '-r300');
% print('-depsc2',sprintf('prof_trap+tri100.eps'), '-r300');

break

% Null
% ****

band = [3.25 4.25];%[3.5 4.1];%
bandtick = [3.25 3.5 3.75 4 4.25];%[3.5 3.6 3.7 3.8 3.9 4.0 4.1];%
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];

figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'YScale','log')
set(gca,'ylim',[.0001 1])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('Null Depth  (log_{10})','FontSize',fsz*1.2)
%y1=plot(lb_t,nul_tmp1,'r:','linewidth',lwz);
%y2=plot(lb_t,nul_tmp2,'g-.','linewidth',lwz);
y3=plot(lb_t,nul_tmp3,'b--','linewidth',lwz);
%y4=plot(lb_t,nul_tmp4,'k','linewidth',lwz);

% %c=plot(lb_t,lb_t.*0+mean(nul_res_sp_b),'k-.','linewidth',lwz)
% %l=legend(' Null Depth (AGPM3)',' Mean null depth');
% legendL1=sprintf(' L1: F=%3.2f h=%3.1fµm slope=%3.2f° R_{filter}=%4.0f',FL1,hL1,penL1,RL1_filter);
% legendL2=sprintf(' L2: F=%3.2f h=%3.1fµm slope=%3.2f° R_{filter}=%4.0f',FL2,hL2,penL2,RL2_filter);
% legendL3=sprintf(' L3: F=%3.2f h=%3.1fµm slope=%3.2f° R_{filter}=%4.0f',FL3,hL3,penL3,RL3_filter);
% legendL4=sprintf(' L4: F=%3.2f h=%3.1fµm slope=%3.2f° R_{filter}=%4.0f',FL4,hL4,penL4,RL4_filter);
% l=legend(legendL1,legendL2,legendL3,legendL4);
% %l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
% set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
% set(l,'box','on','linewidth',2)
% 
% t = title(['     Filter WIDE']);
% set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
% 
% %print('-dpng', [sprintf('rejection_4L_alpha=%3.2f.png',pente_deg)], '-r300');
% print('-dpng', [sprintf('rejection_4L_alpha_variable.png')], '-r300');












