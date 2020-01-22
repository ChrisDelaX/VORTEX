% %Algorythme d'optimisation AGPM
% %------------------------------
% %------------------------------
% %------------------------------
% 
% 
% %Gestion de la mémoire
% %---------------------
% %---------------------
% 
% %Libération
% %----------
% 
% save backup.mat
% clear all;close all;
% warning off MATLAB:singularMatrix
% warning off MATLAB:break_outside_of_loop
% 
% global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
% global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
% 
% %Paramètres de calculs
% %---------------------
% 
% %Troncature en X (N donc 2N+1 ordres au total)
% N=8;%12;%
% 
% %Allocation mémoire
% %------------------
% 
% kx_mat=zeros(2*N+1);
% ky_mat=zeros(2*N+1);
% kIz_mat=zeros(2*N+1);
% kIIIz_mat=zeros(2*N+1);
% E=zeros(2*N+1);
% A=zeros(2*N+1);
% B=zeros(2*N+1);
% D=zeros(2*N+1);
% Omega=zeros(2*(2*N+1));
% Tuu_l=zeros(2*(2*N+1));
% Rud_l=zeros(2*(2*N+1));
% Rdu_l=zeros(2*(2*N+1));
% Tdd_l=zeros(2*(2*N+1));
% F_lp1=zeros(4*(2*N+1));
% X_lp1=zeros(2*(2*N+1));
% 
% 
% %Lectures des entrées
% %--------------------
% %--------------------
% 
% %Choix du profil
% %---------------
% 
% %1: profil rectangulaire simple
% %   -> 1 couche (L=1)
% %2: profil rectangulaire avec couche antireflet dans les creux et sur les bosses
% %   -> 3 couches (L=3)
% %3: profil rectangulaire avec couche antireflet sur les bosses uniquement
% %   -> 2 couches (L=2)
% %4: profil rectangulaire avec couche antireflet continue.5
% %   -> 2 couches (L=2)
% %5: profil rectangulaire simple avec couche d'arret continue
% %   -> 2 couches (L=2)
% %6: profil rectangulaire avec couche d'arret continue et antireflet dans les creux et sur les bosses
% %   -> 4 couches (L=4)
% %7: profil rectangulaire avec couche d'arret continue et antireflet sur les bosses uniquement
% %   -> 3 couches (L=3)
% %8: profil rectangulaire avec couche d'arret continue et antireflet continue
% %   -> 3 couches (L=3)
% %9: profil rectangulaire simple avec couche d'arret discontinue
% %   -> 2 couches (L=2)
% %10: profil rectangulaire avec couche d'arret discontinue et antireflet dans les creux et sur les bosses
% %   -> 4 couches (L=4)
% %11: profil rectangulaire avec couche d'arret discontinue et antireflet sur les bosses uniquement
% %   -> 3 couches (L=3)
% %12: profil rectangulaire avec couche d'arret discontinue et antireflet continue
% %   -> 3 couches (L=3)
% 
% %13: profil trapézoidal simple
% %   -> L couches
% %14: profil trapézoidal avec couche antireflet dans les creux, sur les bosses et adhérent aux parois obliques
% %   -> L+2 couches
% %15: profil trapézoidal avec couche antireflet sur les bosses uniquement
% %   -> L+1 couches
% %16: profil trapézoidal avec couche antireflet continue
% %   -> L+1 couches
% %17: profil trapézoidal simple avec couche d'arret continue
% %   -> L+1 couches
% %18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
% %   -> L+2+1 couches
% %19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
% %   -> L+1+1 couches
% %20: profil trapézoidal avec couche d'arret continue et antireflet continue
% %   -> L+1+1 couches
% %21: profil trapézoidal simple avec couche d'arret discontinue
% %   -> L+1 couches
% %22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
% %   -> L+2+1 couches
% %23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
% %   -> L+1+1 couches
% %24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
% %   -> L+1+1 couches
% %25: profil LETI 1 (sans couche d'arret)
% %26: profil LETI 2 (avec couche d'arret)
% 
% p=13;%1;%3;%
% 
% %Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
% L=25;%50;%8;%
% pente_deg=3.2;%2.75;% 
% pente=rad(pente_deg);
% 
% %pente_min = rad(4);
% %pente_max = rad(7);
% %pente=atan(0.1);%10/180*pi;   % => 10% = 5.71°
% %pente=rad(0.57); % =1%
% %pente=rad(0.86); % =1,5%
% %pente=rad(1.15); % =2%;
% 
% 
% %Choix matériau --> Permittivités
% %--------------------------------
% 
% %1: Vide/Air
% %2: CdTe
% %3: Diamant
% %4: Germanium
% %5: Silicium
% %6: ZnSe
% %7: YF3
% %8: Manuel
% %9: AsGa
% %10: n-laf32
% %11: GASIR 2
% %12: GASIR 1
% %13: InP
% %14: Infrasil
% %15: KRS-5
% %16: Si3N4
% %17: ZnS
% %18: n-lasf44
% %19: ZnSe (T)
% %20: LAH83
% 
% E_man=1.5;
% 
% %Permittivités milieu extérieur incident
% EI_choice=1;
% 
% %Permittivités milieu extérieur émergent
% EIII_choice=3;
% 
% %Permittivités du réseau: EII
% E1_choice=1;  %(la plus basse)
% E2_choice=3;
% 
% %Permittivités de la couche antireflet
% E_AR_choice=1;%14;
% E_AR1_choice=1;
% E_AR2_choice=1;
% 
% %Permittivités de la couche d'arret
% E_arret_choice=1;
% 
% 
% %Paramètres du réseau
% %--------------------
% 
% %Pas du réseau (sublambda!)
% % --> limite Diamant bande N : 3.78  L : 1.4277  K : 0.8391
% Lb=1.42;%4.6;%
% Lb_min=4.7;Lb_max=4.75;
% 
% %Epaisseurs µm
% d=4.2;%4.3;%
% d_min=12;%12.5;%
% d_max=16;%14.5;%
% d_AR=0.336;
% d_AR_min=0.321;d_AR_max=0.351;
% 
% %double AGPM
% %d=d/2;
% %d_min=d_min/2;
% %d_max=d_max/2;
% 
% %Facteurs de remplissage
% F=.36;%0.4;%
% F_min=0.35;%0.42;%
% F_max=0.45;%0.50;%
% 
% %Onde incidente
% %--------------
% 
% %Longueur d'onde µm
% lb=3.5;%11;%      35 9.5 738
% lb_min=3.25;%3.5;%
% lb_max=4.25;%4.1;%
% nlb=81;%11;%
% 
% %Angle incidence non conique
% theta=rad(0);
% theta_min=rad(30); theta_max=rad(50);
% 
% %Angle incidence conique
% phi=rad(0);
% phi_min=rad(0); phi_max=rad(0);
% 
% %Angle polarisation
% psi=rad(45);
% psi_min=rad(0); psi_max=rad(0);
% 
% 
% %Optimisation
% %------------
% 
% 
% % F=0.4;
% % d=4.3;%5.05;
% % x0=[d];
% % [x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
% % x
% 
% % valeurs correctes: 5.1979    0.4651
% % valeurs incorrectes: 4.3031    0.4026
% 
% disp('start')
% datestr(now)
% tic
% 
% x0=[];
% 
% % interpolation linéaire de la pente alpha
% x1= 4.7;
% y1= 2.9;
% x2= 5.57;
% y2= 0.85*rad2deg(atan((1-0.41)*1.42/2/x2));
% 
% 
% pen_a = (y2-y1)/(x2-x1);
% pen_b = y1 - pen_a*x1;
% 
% % Filters
% % -------
% 
% load filters.mat
% x1 = filter_wide(:,1);
% y1 = filter_wide(:,2);
% p_wide = fit(x1,y1,'smoothingspline');
% 
% 
% % Best case (for AGPML4)
% % ----------------------
% 
% F=0.45;%
% d=5.2;%
% %calcul pente
% % pen_interp = d*pen_a + pen_b; 
% % pen_max = rad2deg(atan((1-F)*Lb/2/d));
% % penL4= min(pen_interp,pen_max);
% penL4=2.95;%2.9;%
% pente=rad(penL4);
% % calcul réjection
% [fval1,fval2] = null(x0);
% nul_tmp0 = nul_res_sp_b;
% FL0=F;
% hL0=d;
% %RL0 = 1/mean(nul_tmp4)%*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
% RL0_filter = 1/ (sum(nul_tmp0.*p_wide(lb_t)')/sum(p_wide(lb_t)))
% Trans_tmp0 = (tmp1+tmp2)./2;%.*100
% TransL0 = sum(Trans_tmp0.*p_wide(lb_t)')/sum(p_wide(lb_t))
% TransL0_filter = sum(Trans_tmp0.*p_wide(lb_t)'.*p_wide(lb_t)')/sum(p_wide(lb_t))
% 
% 
% % AGPM-L4
% % -------
% 
% F=0.41;%
% d=4.7;%
% %calcul pente
% % pen_interp = d*pen_a + pen_b; 
% % pen_max = rad2deg(atan((1-F)*Lb/2/d));
% % penL4= min(pen_interp,pen_max);
% penL4=3.1;%2.973;%2.9;%
% pente=rad(penL4);
% % calcul réjection
% [fval1,fval2] = null(x0);
% nul_tmp4 = nul_res_sp_b;
% FL4=F;
% hL4=d;
% %RL4 = 1/mean(nul_tmp4)%*5e-3   % à 2xlambda/D du centre de l'étoile (décroissance naturelle de la PSF)
% RL4_filter = 1/ (sum(nul_tmp4.*p_wide(lb_t)')/sum(p_wide(lb_t)))
% Trans_tmp4 = (tmp1+tmp2)./2;%.*100
% TransL4 = sum(Trans_tmp4.*p_wide(lb_t)')/sum(p_wide(lb_t))
% TransL4_filter = sum(Trans_tmp4.*p_wide(lb_t)'.*p_wide(lb_t)')/sum(p_wide(lb_t))
% 
% TransBS=n2lb.*(2./(n2lb+1)).^2;%.*100    % backside
% 
% 
% disp('end')
% datestr(now)
% toc
% 
% %---------
% 
% 
% %mesures
% %-------
% 
% 
% %transmes = [3.50017500900000,0.901279707000000;3.54987575400000,0.900362319000000;3.59971202300000,0.897018970000000;3.64963503600000,0.889089270000000;3.69959304500000,0.882882883000000;3.74953130900000,0.893693694000000;3.80083618400000,0.865886589000000;3.85059684300000,0.845184518000000;3.90015600600000,0.834532374000000;3.94944707700000,0.833633094000000;4,0.837230216000000;4.05022276200000,0.842625899000000;4.10004100000000,0.843525180000000;];
% transmes = [3.50017500900000,0.901279707000000;3.54987575400000,0.900362319000000;3.59971202300000,0.897018970000000;3.64963503600000,0.889089270000000;3.69959304500000,0.882882883000000;3.74953130900000,0.893693694000000;3.80083618400000,0.865886589000000;3.85059684300000,0.845184518000000;3.90015600600000,0.834532374000000;3.94944707700000,0.833633094000000;4,0.837230216000000;4.05022276200000,0.842625899000000;4.10004100000000,0.843525180000000;];
% 
% absmes0 = [3.500175009,0.98364115;
%     3.69959304500000,0.993247997;
%     3.90015600600000,0.888722592;
%     4.10004100000000,0.874872506;];
% 
% absmes = [3.5,0.985;%0.995;
%     3.6,0.91;%0.92;
%     3.7,1.04;
%     3.8,1.01;
%     3.9,0.76;
%     4.0,0.93;%0.91;
%     4.1,0.93;%0.91;
%     4.2,0.94];%0.92;
% abspol=polyval(polyfit(absmes(:,1),absmes(:,2),3),lb_t);
% 
% mes_ar=[3.50017500900000,0.803817408000000;3.69959304500000,0.816328829000000;3.90015600600000,0.734215470000000;4.10004100000000,0.722866003000000;];
% mes_smooth=[3.50017500900000,0.701642336;3.69959304500000,0.708880309;3.90015600600000,0.630727763;4.10004100000000,0.620466786;];
% 
% 
% arpol=polyval(polyfit(absmes(:,1),absmes(:,2),3),lb_t);
% 
% %AGPM-L4
% trans_AR = 0.981; %1-0.019
% trans_smooth = 0.833121403;
% 
% 
% % figure
% % hold on
% % plot(lb_t,arpol)
% % %plot(absmes(:,1),absmes(:,2),'*b-','linewidth',lwz,'MarkerSize',15);
% % plot(absmes0(:,1),absmes0(:,2),'*r-','linewidth',lwz,'MarkerSize',15);
% % band = [3.5 4.1];%[3.5 4.1];%
% % set(gca,'XLim',band)
% % 
% % break
% 
% %---------
% 
% close all
% warning off MATLAB:polyfit
% 
% fnz = 'Arial'; % fontname
% fsz = 28; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2;  % line width
% 
% 
% band = [3.5 4.1];%[3.5 4.1];%
% %bandtick = [3.25 3.5 3.75 4 4.25];%[3.5 3.6 3.7 3.8 3.9 4.0 4.1];%
% liss_lb_t = [lb_t(1):1e-4:lb_t(end)];
% 
% 
% 
% % TRANS without filter
% %--------------------
% 
% 
% figure('name','trans')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',band)%,'xtick',bandtick)
% %set(gca,'YScale','log')
% set(gca,'ylim',[.6 1.05],'ytick',[.65 .7 .75 .8 .85 .9 .95 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
% ylabel('Transmission','FontSize',fsz*1.2)
% y1=plot(lb_t,Trans_tmp4.*trans_AR,'b--','linewidth',lwz*4);
% y2=plot(lb_t,Trans_tmp4.*trans_AR.*abspol,'color',[0 .5 0],'linewidth',lwz*4);
% y3=plot(transmes(:,1),transmes(:,2),'*-','color',[0 .5 0],'linewidth',lwz,'MarkerSize',15);
% y4=plot(lb_t,Trans_tmp4.*trans_smooth.*abspol,'r-','linewidth',lwz);
% y5=plot(lb_t,abspol.*trans_smooth^2,'k-.','linewidth',lwz*4);
% 
% l=legend(' AGPM + ARG (w/o absorption)',' AGPM + ARG (with absorption)',' Measured',' AGPM w/o ARG (with absorption)',' smooth diamond substrate');
% %l=legend(' Optim1: F=0.403 h=4.303µm',' Optim2: F=0.4 h=4.273µm',' AGPM-L1: F=0.36 h=4.2µm',' AGPM-L2: F=0.36 h=3.8µm');
% set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
% set(l,'box','on','linewidth',2)
% 
% % t = title(['     L-band transmission measurements']);
% % set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
% 
% print('-depsc2',sprintf('trans.eps'), '-r300');
% 
% % break
% 
% 
% 
% % Null + ghost
% %--------------
% 
% trans_ghost = (1-Trans_tmp4).* (1-trans_AR);
% 
% figure('name','lognuldpt')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',band)%,'xtick',bandtick)
% set(gca,'YScale','log')
% set(gca,'ylim',[.0001 .01])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
% ylabel('Null Depth  (log_{10})','FontSize',fsz*1.2)
% y1=plot(lb_t,nul_tmp4+trans_ghost,'color',[0 .5 0],'linewidth',lwz*4);
% y2=plot(lb_t,nul_tmp4,'b--','linewidth',lwz*4);
% y3=plot(lb_t,trans_ghost,'r-','linewidth',lwz);
% y4=plot(lb_t,nul_tmp0,'k-.','linewidth',lwz*4);
% 
% l=legend(' N_{AGPM-L4} = N_{theo} + N_{ghost}',' N_{theo}',' N_{ghost}',' N_{theo} (optimal parameters)');
% set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
% set(l,'box','on','linewidth',2)
% 
% % t = title(['     L-band expected performance']);
% % set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
% 
% print('-depsc2',sprintf('null.eps'), '-r300');
% 

% % Profile d'attenuation
% % ---------------------
% 
% close all
% warning off MATLAB:polyfit
% 
% fnz = 'Arial'; % fontname
% fsz = 28; % fontsize
% fwz = 'normal';%'Bold'; % fontweight
% msz = 8; % marker size
% lwz = 2;  % line width
% 
% load labLESIA.mat
% 
% figure('name','lognuldpt')
% set(gcf,'color',[1 1 1])
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
% set(gcf,'Position',[400    0   800   600],'PaperUnits','points','PaperPosition',[400    0   800   600]); % paper 
% hold on
% grid on
% set(gca,'box','on','linewidth',2)
% set(gca,'XLim',[0 3.5],'xtick',[0 1 2 3])
% set(gca,'YScale','log')
% set(gca,'ylim',[.001 1],'ytick',[.01 .1 .5 1])
% set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% set(gca,'XMinorTick','on')
% set(gca,'YMinorTick','on')
% xlabel('Angular separation (\lambda/D)','FontSize',fsz*1.2)
% ylabel('Off-axis transmission','FontSize',fsz*1.2)
% %ylabel('Relative intensity','FontSize',fsz*1.2)
% 
% 
% x = 0:.01:3.5;
% BES_theo = BESSELJ(1,x.*pi./sqrt(3));
% y_theo=1-((4.*BES_theo.^2)./(x.*pi/sqrt(3)).^2);
% y1=plot(x,y_theo,'b','linewidth',lwz*2);
% 
% y2=plot(att_mes(1,1),att_mes(1,2),'s-','color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2);
% 
% xx=att_mes(:,1);
% yy=att_mes(:,2);
% xx(end+1)=xx(end);
% xx(end-1)=xx(end-2)+(xx(end)-xx(end-2))*.45;
% yy(end+1)=yy(end);
% %yy(end-1)=yy(end-2)+(yy(end)-yy(end-2))*.99;
% 
% y_mes=spline(xx,yy,x);
% y3=plot(att_mes(:,1),att_mes(:,2),'s',x,y_mes);
% set(y3,'color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2)
% 
% l=legend(' Theoretical',' Measured');
% set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','best')
% set(l,'box','on','linewidth',2)
% 
% print('-depsc2',sprintf('offaxis.eps'), '-r300');
% %print('-depsc2',sprintf('test_rel_int.eps'), '-r300');



% Profile CORONO
% --------------

close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width

figure('name','lognuldpt')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
grid on
set(gca,'box','on','linewidth',2)
set(gca,'XLim',[0 4.6]);%,'xtick',[0 1 2 3])
set(gca,'YScale','log')
set(gca,'ylim',[.000001 1]);%,'ytick',[1 .01 .0001 .000001])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Angular separation (\lambda/D)','FontSize',fsz*1.2)
ylabel('Relative intensity','FontSize',fsz*1.2)

load coro_profile.mat

offset = 5;
niv_psf = x450PSF_profile0x2B00x2E86(1)+offset;
xx = 0:.1:5;
xx=xx./.8;
ypsf =  (x450PSF_profile0x2B00x2E86+offset)./niv_psf;
ycoro_mean = x450Mean_coro_profile0x2B00x2E86./niv_psf;
ycoro_std = x450stddev_coro_profile0x2B00x2E86./niv_psf;
ybkg = Level./niv_psf;

y1=plot(xx,ypsf,'b','linewidth',lwz*2);
y2=plot(xx,ycoro_mean,'--','color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2);

clear zone_std
zone_std(:,1)=(ycoro_mean-ycoro_std);;
zone_std(find(zone_std<0))=1e-6;
zone_std(:,2)=(ycoro_mean+ycoro_std-zone_std(:,1))

h=area(xx',zone_std)
set(h(1),'visible','off')
set(h,'LineStyle',':','EdgeColor',[0 .5 0],'linewidth',lwz*2)

y2=plot(xx,ycoro_mean,'--','color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2);

%y4=plot(xx,ycoro_mean+ycoro_std,':','color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2);
%y5=plot(xx,ycoro_mean-ycoro_std,':','color',[0 .5 0],'MarkerSize',15,'linewidth',lwz*2);


l=legend(' Off-axis profile',' Coronagraphic profile');
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','southwest')
set(l,'box','on','linewidth',2)



print('-depsc2',sprintf('psf_coro.eps'), '-r300');






