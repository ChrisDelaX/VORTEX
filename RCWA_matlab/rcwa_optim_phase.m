%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Gestion de la mémoire
%---------------------
%---------------------

%Libération 
%----------

%save data/backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global ordre

%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)  
N=12;%6;%

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

p=13;

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=160;
pente=rad(5);
pente_max=rad(5.9);
pente_min=rad(0.9);

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

E_man=1.494^2;%1.44;

%Permittivités milieu extérieur incident   
EI_choice=1;

%Permittivités milieu extérieur émergent 
EIII_choice=3;%8;%

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;%8;%

%Permittivités de la couche antireflet
E_AR_choice=14;
E_AR1_choice=1;
E_AR2_choice=1;

%Permittivités de la couche d'arret
E_arret_choice=1;


%Paramètres du réseau
%--------------------    

%Pas du réseau (sublambda!)     --> limite Diamant bande N : 3.78
Lb=4.679;%20;%3.6807;%3.505;%
Lb_min=3.957;Lb_max=3.973;

%Epaisseurs
d=15.765;%12;%1.1;%
d_min=13.76;d_max=13.92;
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage         
F=0.2997;%0.72;%0.6025;%
F_min=0.697;F_max=0.701;

%Onde incidente 
%--------------

%Longueur d'onde
lb=.63282;%10;%
lb_min=10.5;%10.5;%
lb_max=10.5;%10.51;%
nlb=1;%11;%

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

% x0=[];
%fld1=(1-F)*Lb;
x0=[];
%x0=[Lb,d];


tic

    ordre = 0;


    npts=120;%3;%
    for i=1:npts

        %Lb=Lb_min+(Lb_max-Lb_min)*(i-1)/(npts-1);
        %d=d_min+(d_max-d_min)*(i-1)/(npts-1);
        %F=F_min+(F_max-F_min)*(i-1)/(npts-1);
        pente=pente_min+(pente_max-pente_min)*(i-1)/(npts-1);

        i
        %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
        %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',35,'Display','iter','TolX',5e-3,'TolFun',1e-5));
        fval=null(x0);


        Xparam(i)=pente/pi*180;
        Yparam(i)=tmp3(1);
        nulldpt(i)=fval;
    end

toc
save rcwa_pente_phase.mat


%MOY TE-TM
%---------
% TEparam=tmp4;
% TMparam=tmp5;
% TEMmoy=mean(tmp4+tmp5)/2;


%MIN
%---
% p = polyfit(Xparam,nulldpt,2);
% nbp = size(p,2);
% for i=1:nbp-1
%     q(1,i)=p(1,i)*(nbp-i);
% end
% Froot = roots(q)


%graphiques :
%------------

%Sans balayage
%-------------
%sansBAL

%Pente
%-----------------------
figure('FileName','Diamant')
hold on
%title('AGPM Diamond/SiO2 : K-band')
%plot(Xparam,YparamTE,'-',Xparam,YparamTM,'--','Linewidth',2)

plot(Xparam,Yparam,'-','Linewidth',2)
xlabel('Pente (°)')
%ylabel('DET_s  et  DET_p')
%legend('TE','TM')
%legend


%DE(N)
%-----
% figure('FileName','Diamant')
% hold on
% %nulldepth = subplot(1,3,1);
% %title('AGPM Diamond/SiO2 : K-band')
% %tmp12=mean((tmp1+tmp2)./2);
% %TsansZOG=2./(n2lb+1);
% plot(lb_t,tmp1./2,'-',lb_t,tmp2./2,'--','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('DET_s  et  DET_p')
% legend('TE','TM')


%Filling Factor
%--------------
% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(1,3,1);
% %title('AGPM Diamond/SiO2 : K-band')
% plot(Xparam,nulldpt,'Marker','.')
% xlabel('Filling factor')
% ylabel('Mean Null Depth')
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


% %Total Grating thickness
% %-----------------------
% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(1,3,1);
% %title('AGPM Diamond/SiO2 : K-band')
% plot(Xparam,nulldpt,'Marker','o')
% xlabel('Total Grating thickness (\mum)')
% ylabel('Mean Null Depth')
% axis square
% set(gca,'YScale','log')
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
