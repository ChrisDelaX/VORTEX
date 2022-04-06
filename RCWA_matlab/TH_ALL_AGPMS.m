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





% % Calcul rcwa
% % -----------
% F=0.45;%0.40;%
% pente=rad(3);
% lb_min=3.5;%11;%3.5;%2;%1.5;%1.15;%.65;%
% lb_max=4.1;%13.2;%4.1;%2.4;%1.8;%1.4;%.78;%
% Lb=1.42;%4.6;0.82;%0.61;%0.48;%0.26;%
% d=5;%14.3957;%2.5830;%1.9215;%1.5022;%3*Lb;%0.8083;%
% nlb=81;
% 
% nlb=11;%
% x0=[d];
% [x,fval,exitflag] = fminsearchbnd(@null,x0,d*.5,d*2,optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
% %x0=[F d];
% %[x,fval,exitflag] = fminsearchbnd(@null,x0,[F/2 d/2],[F*2 d*2],optimset('MaxIter',35,'Display','iter','TolX',1e-5,'TolFun',1e-6));
% x
% [nuldpt,sigerr]=null(x)
% R=1/nuldpt
% break
% 
% F=0.4;
% d=0.8083;%1.5022;%1.9215;%2.5830;%3*Lb;%
% [nuldpt1,sigerr1]=null()
% R1=1/nuldpt1
% nul_tmp1 = nul_res_sp_b;
% F1=F;
% h1=d;
% F=0.4623;%0.4484;%0.4532;%0.4506;    
% d=0.9762;%1.7453;%2.2625;%3.0166;
% [nuldpt2,sigerr2]=null()
% R2=1/nuldpt2
% nul_tmp2 = nul_res_sp_b;
% F2=F;
% h2=d;
% %break
% 
% % Figures
% % -------
% close all
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width 
% bandx=[lb_min lb_max];
% figure('name','null')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800 800],'PaperUnits','points','PaperPosition',[400    0   800 800]); % paper 
% % null depth
% % -----------
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandx)%,'xtick',bandtick)
% set(gca,'YScale','log')
% set(gca,'ylim',[.0001 .01])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength (µm)','FontSize',fsz*1.2)
% ylabel('Null Depth','FontSize',fsz*1.2)
% y1=plot(lb_t,nul_tmp1,'ko-','linewidth',2*lwz);
% y2=plot(lb_t,lb_t.*0+nuldpt1,'k--','linewidth',2*lwz);
% y3=plot(lb_t,nul_tmp2,'o-','color',mygreen,'linewidth',2*lwz);
% y4=plot(lb_t,lb_t.*0+nuldpt2,'--','color',mygreen,'linewidth',2*lwz);
% %leg=legend(' relaxed: F=40%  h=2.58µm  R=1005',' mean null (relaxed)','  optim. : F=45%  h=3.02µm  R=2366',' mean null (optim)');
% %leg=legend(' relaxed: F=40%  h=1.92µm  R=933',' mean null (relaxed)','  optim. : F=45%  h=2.26µm  R=2341',' mean null (optim)');
% %leg=legend(' relaxed: F=40%  h=1.50µm  R=962',' mean null (relaxed)','  optim. : F=45%  h=1.75µm  R=2353',' mean null (optim)');
% leg=legend(' relaxed: F=40%  h=0.81µm  R=732',' mean null (relaxed)','  optim. : F=45%  h=0.98µm  R=2144',' mean null (optim)');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% %tit=title(' Ks-band [2-2.4 µm]: \alpha = 3.0° and \Lambda = 0.82 µm');
% %tit=title(' H-band [1.5-1.8 µm]: \alpha = 3.0° and \Lambda = 0.61 µm');
% %tit=title(' J-band [1.15-1.4 µm]: \alpha = 3.0° and \Lambda = 0.48 µm');
% tit=title(' visible [0.65-0.78 µm]: \alpha = 3.0° and \Lambda = 0.26 µm');
% set(tit,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')%,'horizontalalignment','left')
% %print('-depsc2',sprintf('null_opt_Ks.eps'), '-r300');
% %print('-depsc2',sprintf('null_opt_H.eps'), '-r300');
% %print('-depsc2',sprintf('null_opt_J.eps'), '-r300');
% print('-depsc2',sprintf('null_opt_vis.eps'), '-r300');



% break
 
% % Figures
% % -------
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 28; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width 
% bandx=[lb_min lb_max];
% close all
% figure('name','null')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   1600 800],'PaperUnits','points','PaperPosition',[400    0   2000 800]); % paper 
% % trans
% % ----------
% fig3 = subplot(2,3,1);
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[0.8 1])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz*1.2)
% ylabel('Transmittance','FontSize',fsz*1.2)
% plot(lb_t,tmp1,'--','color',myblue,'linewidth',lwz)
% plot(lb_t,tmp2,'--','color',myred,'linewidth',lwz)
% plot(lb_t,(tmp1+tmp2)./2,'k-','linewidth',lwz)
% %plot(lb_t,(n2lb-1)./(n2lb+1),'k.:','linewidth',lwz)
% leg=legend(' T_{TE}(\lambda)',' T_{TM}(\lambda)',' mean T');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','south')
% set(leg,'box','on','linewidth',lwz)
% % refl
% % -----------
% fig4 = subplot(2,3,4);
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandx)%,'xtick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16])
% %set(gca,'YScale','log')%,'XMinorGrid','on','XAxisLocation','top')
% set(gca,'ylim',[0 0.2])%,'ytick',[1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1])%,'ydir','reverse')
% set(gca,'Fontname',fnz,'FontSize',0.9*fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength   \lambda  (µm)','FontSize',fsz*1.2)
% ylabel('Reflectance','FontSize',fsz*1.2)
% plot(lb_t,tmp4,'--','color',myblue,'linewidth',lwz)
% plot(lb_t,tmp5,'--','color',myred,'linewidth',lwz)
% plot(lb_t,(tmp4+tmp5)./2,'k-','linewidth',lwz)
% leg=legend(' R_{TE}(\lambda)',' R_{TM}(\lambda)',' mean R');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% % null depth
% % -----------
% fig4 = subplot(2,3,[2:3 5:6]);
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',bandx)%,'xtick',bandtick)
% set(gca,'YScale','log')
% set(gca,'ylim',[.0001 .1])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength (µm)','FontSize',fsz*1.2)
% ylabel('Null Depth','FontSize',fsz*1.2)
% y1=plot(lb_t,nul_res_sp_b,'o-','color',mygreen,'linewidth',2*lwz);
% y2=plot(lb_t,nuldpt,'k--','linewidth',2*lwz);
% legendK1=sprintf(' Ks-band AGPM: F=%3.2f h=%3.1fµm R=%4.0f',F,d,1./nuldpt);
% leg=legend(legendK1,' mean null depth');
% %l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('null_Ks.eps'), '-r300');


% break

% %Optimisation
% %------------
% Lb=4.6;
% lb_min=11;%
% lb_max=13.2;%
% nlb=81;%11;%
% disp('start')
% datestr(now)
% tic
% x0=[];
% % AGPM-N1
% % -------
% F=0.3152;%0.35;%0.36;%
% d=14.75;%4.2;%4.5;%
% penN1=2.75;%2.4;%2.8;%2.7;%2.9;%
% pente=rad(penN1);
% % calcul réjection
% [fval1,fval2] = null(x0);
% nul_tmp1 = nul_res_sp_b;
% FN1=F;
% hN1=d;
% % AGPM-N3
% % -------
% F=0.3587;%0.44;%
% d=12.6;%5.8;%6;%6;%
% penN3=2.75;%3.25;%2.9;%pen_max;%
% pente=rad(penN3);
% % calcul réjection
% [fval1,fval2] = null(x0);
% nul_tmp3 = nul_res_sp_b;
% FN3=F;
% hN3=d;
% % AGPM-N4
% % -------
% F=0.3967;%0.41;%
% d=13.8;%4.7;%
% penN4=2.75;%3.0;%2.9;%
% pente=rad(penN4);
% % calcul réjection
% [fval1,fval2] = null(x0);
% nul_tmp4 = nul_res_sp_b;
% FN4=F;
% hN4=d;
% disp('end')
% datestr(now)
% toc
% 
% 
% % Figures
% % -------
% close all
% warning off MATLAB:polyfit
% mygreen = [0 .5 0];
% myred = [1 .2 0];
% myblue = [0 .2 1];
% fnz = 'Arial'; % fontname
% fsz = 26; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2.2;  % line width 
% 
% band = [11 13.2];%[3.25 4.25];%
% %bandtick = [3.25 3.5 3.75 4 4.25];%[3.5 3.6 3.7 3.8 3.9 4.0 4.1];%
% liss_lb_t = [lb_t(1):1e-4:lb_t(end)];
% figure('name','lognuldpt')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',band)%,'xtick',bandtick)
% set(gca,'YScale','log')
% set(gca,'ylim',[.0001 1])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength (µm)','FontSize',fsz*1.2)
% ylabel('Null Depth','FontSize',fsz*1.2)
% y1=plot(lb_t,nul_tmp1,'--','color',myred,'linewidth',lwz);
% y3=plot(lb_t,nul_tmp3,'k-','linewidth',2*lwz);
% y4=plot(lb_t,nul_tmp4,'-','color',mygreen,'linewidth',2*lwz);
% legendN1=sprintf(' N1: F=%3.2f h=%3.1fµm R=%4.0f',FN1,hN1,1./mean(nul_tmp1));
% legendN3=sprintf(' N3: F=%3.2f h=%3.1fµm R=%4.0f',FN3,hN3,1./mean(nul_tmp3));
% legendN4=sprintf(' N4: F=%3.2f h=%3.1fµm R=%4.0f',FN4,hN4,1./mean(nul_tmp4));
% leg=legend(legendN1,legendN3,legendN4);
% %l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
% set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
% set(leg,'box','on','linewidth',lwz)
% print('-depsc2',sprintf('null_N_ALL.eps'), '-r300');


%break
%break


%Optimisation
%------------
Lb=1.42;
lb_min=3.5;%
lb_max=4.1;%
nlb=81;%11;%
disp('start')
datestr(now)
tic
x0=[];
% AGPM-L1
% -------
F=0.35;%0.36;%
d=4.2;%4.5;%
penL1=2.4;%2.8;%2.7;%2.9;%
pente=rad(penL1);
% calcul réjection
[fval1,fval2] = null(x0);
nul_tmp1 = nul_res_sp_b;
FL1=F;
hL1=d;
% AGPM-L2
% -------
F=0.36;%
d=3.6;%
penL2=2.65;%3;%2.9;%
pente=rad(penL2);
% calcul réjection
[fval1,fval2] = null(x0);
nul_tmp2 = nul_res_sp_b;
FL2=F;
hL2=d;
% AGPM-L3
% -------
F=0.44;%
d=5.8;%6;%6;%
penL3=3.25;%2.9;%pen_max;%
pente=rad(penL3);
% calcul réjection
[fval1,fval2] = null(x0);
nul_tmp3 = nul_res_sp_b;
FL3=F;
hL3=d;
% AGPM-L4
% -------
F=0.41;%
d=4.7;%
penL4=3.0;%2.9;%
pente=rad(penL4);
% calcul réjection
[fval1,fval2] = null(x0);
nul_tmp4 = nul_res_sp_b;
FL4=F;
hL4=d;
disp('end')
datestr(now)
toc


% Figures
% -------
close all
warning off MATLAB:polyfit
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
fnz = 'Arial'; % fontname
fsz = 26; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2.2;  % line width 

band = [3.5 4.1];%[3.25 4.25];%
%bandtick = [3.25 3.5 3.75 4 4.25];%[3.5 3.6 3.7 3.8 3.9 4.0 4.1];%
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];
figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'box','on','linewidth',2)
set(gca,'XLim',band)%,'xtick',bandtick)
set(gca,'YScale','log')
set(gca,'ylim',[.0001 1])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength (µm)','FontSize',fsz*1.2)
ylabel('Null Depth','FontSize',fsz*1.2)
y1=plot(lb_t,nul_tmp1,'--','color',myred,'linewidth',lwz);
y2=plot(lb_t,nul_tmp2,'--','color',myblue,'linewidth',lwz);
y3=plot(lb_t,nul_tmp3,'k-','linewidth',2*lwz);
y4=plot(lb_t,nul_tmp4,'-','color',mygreen,'linewidth',2*lwz);
legendL1=sprintf(' L1: F=%3.2f h=%3.1fµm angle=%3.2f° R=%4.0f',FL1,hL1,penL1,1./mean(nul_tmp1));
legendL2=sprintf(' L2: F=%3.2f h=%3.1fµm angle=%3.2f° R=%4.0f',FL2,hL2,penL2,1./mean(nul_tmp2));
legendL3=sprintf(' L3: F=%3.2f h=%3.1fµm angle=%3.2f° R=%4.0f',FL3,hL3,penL3,1./mean(nul_tmp3));
legendL4=sprintf(' L4: F=%3.2f h=%3.1fµm angle=%3.2f° R=%4.0f',FL4,hL4,penL4,1./mean(nul_tmp4));
leg=legend(legendL1,legendL2,legendL3,legendL4);
%l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
set(leg,'Fontname',fnz,'FontWeight',fwz,'FontSize',fsz,'location','north')
set(leg,'box','on','linewidth',lwz)
print('-depsc2',sprintf('null_L_ALL.eps'), '-r300');












