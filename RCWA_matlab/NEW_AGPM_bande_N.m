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

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 limZOG

%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)  
N=12;%8;%20;%

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

p=30;%3;%

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=15;%30;%
pente=rad(4.3);
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

E_man=1.494;%1.44;

%Permittivités milieu extérieur incident   
EI_choice=1;

%Permittivités milieu extérieur émergent 
EIII_choice=3;

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;

%Permittivités de la couche antireflet
E_AR_choice=14;
E_AR1_choice=1;
E_AR2_choice=1;

%Permittivités de la couche d'arret
E_arret_choice=1;


%Paramètres du réseau
%--------------------    

%Pas du réseau (sublambda!)     
% --> limite Diamant bande N : 3.78  L : 1.4697  K : 0.8391
Lb=1.47;
Lb_min=4.7;Lb_max=4.75;

%Epaisseurs
d=9;%8.87;%
d_min=2;d_max=7;
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage         
F=0.3547;
F_min=0.3;F_max=0.7;

%Onde incidente 
%--------------

%Longueur d'onde
lb=11;%0.55;%      35 9.5 738
lb_min=3.5;%0.3;%
lb_max=4.1;
nlb=11;%41;%

%Angle incidence non conique
theta=rad(0);
theta_min=rad(30); theta_max=rad(50);

%Angle incidence conique
phi=rad(0);
phi_min=rad(0); phi_max=rad(0);

%Angle polarisation
psi=rad(45);
psi_min=rad(0); psi_max=rad(0);


% Optim Anti-reflet :
% -------------------
% m = 7; %entier impair
% lb = lb_min;
% material
% n1 = sqrt(EI);
% n3 = sqrt(EIII);
% n2 = sqrt(n1*n3);
% d_min=m*lb/4/n2
% [x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
% lb = lb_max;
% material
% n1 = sqrt(EI);
% n3 = sqrt(EIII);
% n2 = sqrt(n1*n3);
% d_max=m*lb/4/n2


%Optimisation
%------------

 x0=[];
%fld1=(1-F)*Lb;
%x0=[F, d_AR];

tic
npts=11;%1;%
for i=1:npts


    %Lb=Lb_min+(Lb_max-Lb_min)*(i-1)/(npts-1);
    d=d_min+(d_max-d_min)*(i-1)/(npts-1);
    %F=F_min+(F_max-F_min)*(i-1)/(npts-1);

    i
    %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-6,'TolFun',1e-6));
    %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',35,'Display','iter','TolX',1e-6,'TolFun',1e-6))
    %[x,fval,exitflag] = fminsearch(@null_opt,x0,optimset('MaxIter',35,'Display','iter','TolX',1e-6,'TolFun',1e-6));

    fval=null([]);
    

    
    Xparam(i)=d;
    %Yparam(i,:)=x;
    nulldpt(i)=null_res_sp;
end
toc
save NEW_AGPM_bande_N.mat

%save AGPM_Lband_alpha=0_d=4.7.mat
% 
% TEparam=tmp4;
% TMparam=tmp5;
% TEMmoy=mean(tmp4+tmp5)/2;


%MIN
%---
% p = polyfit(Xparam,nulldpt,4);
% nbp = size(p,2);
% for i=1:nbp-1
%     q(1,i)=p(1,i)*(nbp-i);
% end
% Froot = roots(q)


%graphiques :
%------------

%Filling Factor
%--------------

%Yparam(:,2)=Yparam(:,1);
%Yparam(:,1)= Lb;

% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(1,3,1);
% %title('AGPM Diamond/SiO2 : K-band')
% plot(Xparam,nulldpt,'Marker','.')
% xlabel('Depth (\mum)')
% title('Null Depth @ 3.8 \mum  |  Profile slope \alpha = 4.3 °')
% axis square
% set(gca,'YScale','log')
% % ParamY1 = subplot(1,4,2);
% % plot(Xparam,Yparam(:,1),'Marker','.')
% % xlabel('Depth (\mum)')
% % ylabel('Period (\mum)')
% ParamY2 = subplot(1,3,2);
% plot(Xparam,Yparam(:,1),'Marker','.')
% xlabel('Depth (\mum)')
% title('Filling factor')
% AspRatio = subplot(1,3,3);
% for i=1:npts
%     if Yparam(i) < 0.5
%         F_AR(i) = Yparam(i);
%     else
%         F_AR(i) = 1-Yparam(i);
%     end
% end
% AR = Xparam./(F_AR.*Lb);
% plot(Xparam,AR,'Marker','.')
% xlabel('Depth (\mum)')
% title('Aspect Ratio')


% %Total Grating thickness
% %-----------------------
figure('FileName','Diamant')
hold on
nulldepth = subplot(1,3,1);
title('AGPM Diamond proto N-band')
plot(Xparam,nulldpt,'Marker','o')
xlabel('Total Grating thickness (\mum)')
ylabel('Mean Null Depth')
axis square
set(gca,'YScale','log')
% ParamY1 = subplot(1,3,2);
% plot(Xparam,Yparam(:,1),'Marker','o')
% xlabel('Total Grating thickness (\mum)')
% ylabel('Period (\mum)')
% ParamY2 = subplot(1,3,3);
% plot(Xparam,Yparam(:,2),'Marker','.')
% xlabel('Total Grating thickness (\mum)')
% ylabel('Filling factor')



% %Period
% %-----------------------
% figure('FileName','Diamant')
% hold on
% %title('AGPM Diamond/SiO2 : K-band')
% plot(Xparam,nulldpt,'Marker','o')
% xlabel('Period (\mum)')
% ylabel('Mean Null Depth')
% axis square
% set(gca,'YScale','log')



%Sans balayage
%-------------

%sansbal

% figure('FileName','Diamant')
% hold on
% %grid on
% nulldepth = subplot(2,2,1);
% plot(lb_t,nul_res_sp_b,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
% axis square
% set(nulldepth,'YScale','log')
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Null Depth')
% title('AGPM Diamond : N-band')
% gtext({['Slope \alpha: ',num2str(pente/pi*180),' °'],['Period \Lambda: ',num2str(Lb),' \mum'],['Total Grating thickness d: ',num2str(d),' \mum'],['Filling factor F: ',num2str(F*100),' %']})
% transmittance = subplot(2,2,3);
% tmp12=mean(tmp1+tmp2)/2;
% plot(lb_t,tmp1,'-',lb_t,tmp2,'--',lb_t,tmp12,'k:','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Transmittance')
% legend('TE','TM')
% phaseshift = subplot(2,2,2);
% tmp3moy=mean(tmp3);
% plot(lb_t,tmp3,lb_t,tmp3moy,'k:',lb_t,pi,'k-.','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('TE-TM Phase Shift  \Delta\Phi_{TE-TM}  (rad)')
% reflectance = subplot(2,2,4);
% tmp45=mean(tmp4+tmp5)/2;
% plot(lb_t,tmp4,'-',lb_t,tmp5,'--',lb_t,tmp45,'k:','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Reflectance')
% legend('TE','TM')


%MEAN
%----
% figure 
% p = polyfit(lb_t,nul_res_sp_b,16);
% nptgr=100;
% fin=size(lb_t,2);
% ecart=(lb_t(fin)-lb_t(1))/(nptgr-1);
% for i=1:nptgr
%     x(i)=lb_t(1)+ecart*(i-1);
% end
% y=polyval(p,x);
% null_res_sp_2=mean(y)
% plot(x,y,'k.-',lb_t,null_res_sp_2,'k--','Linewidth',2)
% set(gca,'YScale','log')



% longueur d'onde
% ---------------

% 
% figure('FileName','n3')
% hold on
% plot(lb_t,n3,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
% axis square
% xlabel('Wavelength \lambda (\mum)')
% ylabel('n3')
% 
% figure('FileName','n2')
% hold on
% plot(lb_t,n2,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
% axis square
% xlabel('Wavelength \lambda (\mum)')
% ylabel('n2')
% 
% 
% figure('FileName','Fnew')
% hold on
% plot(lb_t,Fnew,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
% axis square
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Filling Factor')
% 
% figure('FileName','dnew')
% hold on
% plot(lb_t,dnew,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
% axis square
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Depth')


% periode
% -------
% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(1,3,1);
% tmp45=mean(tmp4+tmp5)/2;
% %plot(Xparam,TEparam,'-',Xparam,TMparam,'--',Xparam,TEMmoy,'k:','Linewidth',2)
% 
% xlabel('Period \Lambda (\mum)')
% ylabel('Reflectance')
% axis square


% set(gca,'YScale','log')
% ParamY1 = subplot(1,3,2);
% plot(Xparam,Yparam(:,1),'Marker','.')
% xlabel('Filling factor')
% ylabel('Period (\mum)')
% ParamY2 = subplot(1,3,3);
% plot(Xparam,Yparam(:,2),'Marker','.')
% xlabel('Filling factor')
% ylabel('Total Grating thickness (\mum)')
