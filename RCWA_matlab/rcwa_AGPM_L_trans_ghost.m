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

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof ghost lbtmp

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
L=25;%50;%8;%
pente_deg=2.75; %4.3 ou 5.8
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
EI_choice=1;%3;%

%Permittivités milieu extérieur émergent
EIII_choice=3;%1;%

%Permittivités du réseau: EII
E1_choice=1;%3;%  %(la plus basse)
E2_choice=3;%1;%

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
lb_min=3.5;%3.3;%11;%
lb_max=4.1861386;%4.3;%4.1;%5;%13.2;%
nlb=90;%100;%150;%11;%61;%81;%

lbtmpfct
nlb = size(lbtmp,1);

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

% 
% % optim1
% F=0.47;%0.4026;%
% d=5.2;%4.3031;%
% [fval1,fval2] = null(x0);
% mean(nul_res_sp_b)*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
% nul_opt1 = nul_res_sp_b;
% 
% % optim2
% F=0.4;
% d=4.3;%4.2731;%
% [fval1,fval2] = null(x0);
% mean(nul_res_sp_b)*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
% nul_opt2 = nul_res_sp_b;




% AGPM-L1
F=0.36;
d=4.2;
[fval1,fval2] = null(x0);
mean(nul_res_sp_b)%*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
dphi_tmp1 = tmp3;
transL1=(tmp1+tmp2)./2;%.*100
mattemp(:,5) = nul_res_sp_b;


% AGPM-L2
F=0.36;
d=3.8;%
[fval1,fval2] = null(x0);
mean(nul_res_sp_b)%*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
dphi_tmp2 = tmp3;
transL2=(tmp1+tmp2)./2;%.*100

transBS=n2lb.*(2./(n2lb+1)).^2;%.*100    % backside


mattemp(:,1) = lb_t .* 1000;
mattemp(:,2) = transL1';
mattemp(:,3) = transL2';
mattemp(:,4) = transBS';
mattemp(:,6) = nul_res_sp_b;



%break 


disp('end')
datestr(now)
toc

%---------



close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width


band = [3.5 4.1];%[11 13.2];%
bandtick = [3.5 3.6 3.7 3.8 3.9 4.0 4.1];%[11.0 11.5 12.0 12.5 13.0];%
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];



figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',band,'xtick',bandtick)
%set(gca,'YScale','log')
%set(gca,'ylim',[.0001 .2])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('Transmission (%)','FontSize',fsz*1.2)
y1=plot(lb_t,transL1,'b','linewidth',lwz);
y2=plot(lb_t,transL2,'b:','linewidth',lwz);
y3=plot(lb_t,transBS,'r','linewidth',lwz);

%c=plot(lb_t,lb_t.*0+mean(nul_res_sp_b),'k-.','linewidth',lwz)
%l=legend(' Null Depth (AGPM3)',' Mean null depth');

l=legend(' AGPM-L1',' AGPM-L2',' Backside (without AR)');
%l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)

t=title(['     NACO AGPMs (L-band)']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')


%print('-dpng', ['lognuldpt_L_trans.png'], '-r300');


% l=legend(' AGPM-L1',' AGPM-L2',' without SWG');
% set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
% set(l,'box','on','linewidth',2)
% 
% print('-dpng', ['test.png'], '-r300');
















