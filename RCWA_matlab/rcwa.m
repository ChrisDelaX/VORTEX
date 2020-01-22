%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Gestion de la mémoire
%---------------------
%---------------------

%Libération 
%----------

save data/backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 dphi_sp_T nul_res_sp_b null_res_sp

%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)  
N=8;

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

p=3;

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=1;


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

E_man=1.44;
%LaF2: E_man=1.57;
%BaF2: E_man(l=1.97 microns)=1.4647
%BaF2: E_man(l=2.3253 microns)=1.4635
%BaF2: E_man(l=2.6738 microns)=1.4623

%Permittivités milieu extérieur incident   
EI_choice=1;

%Permittivités milieu extérieur émergent 
EIII_choice=3;

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;

%Permittivités de la couche antireflet
E_AR_choice=14;
E_AR1_choice=1; %7;
E_AR2_choice=1; %17;

%Permittivités de la couche d'arret
E_arret_choice=1; %16;


%Paramètres du réseau
%--------------------    

%Pas du réseau (sublambda!)
Lb=0.768;%0.737;%0.59;%1.5;         %bande K (2-2,4) SiO2 ZOG -> Lb_max (1,3906-1,6762)
Lb_min=0.4;Lb_max=0.8;

%Epaisseurs
d=3.223;%1.11;%1.13;%1.11;%1.2;
d_min=3.19;d_max=3.22;
d_AR=0.336;%0.224;%0.25;%0.35;
d_AR_min=0;
d_AR_max=0.6;
d_AR1=0.5; d_AR2=0.5;
d_arret=0.1;


%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage         
F=0.663;%0.68;%0.55;%0.531;    % 0.85/1.6
F_min=0.65;F_max=0.7;

%Onde incidente 
%--------------

%Longueur d'onde
lb=2.2;
lb_min=2;lb_max=2.4;
nlb=21;

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

 x0=[];
%x0=[Lb,d,d_AR];
% x0=[Lb,d,d_AR,F];

tic
npts=1;
for i=1:npts

    %Lb=Lb_min+(Lb_max-Lb_min)*(i-1)/(npts-1);
    %d=d_min+(d_max-d_min)*(i-1)/(npts-1);
    %d_AR=d_AR_min+(d_AR_max-d_AR_min)*(i-1)/(npts-1);
    %F=F_min+(F_max-F_min)*(i-1)/(npts-1);

    i
    %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
    %[x,fval,exitflag] = fminsearch(@null,x0,optimset('MaxIter',30,'Display','iter','TolX',5e-3,'TolFun',1e-5));
    fval=null(x0);

    
    %Xparam(i)=F;
    %Yparam(i,:)=x;
    nulldpt(i)=fval;
end
toc


%graphiques :
%------------

%Filling Factor
%--------------
% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(2,2,1);
% title('AGPM Diamond/SiO2 : K-band')
% plot(Xparam,nulldpt,'Marker','.')
% xlabel('Filling factor')
% ylabel('Mean Null Depth')
% axis square
% set(gca,'YScale','log')
% ParamY1 = subplot(2,2,2);
% plot(Xparam,Yparam(:,1),'Marker','.')
% xlabel('Filling factor')
% ylabel('Period (\mum)')
% ParamY2 = subplot(2,2,3);
% plot(Xparam,Yparam(:,2),'Marker','.')
% xlabel('Filling factor')
% ylabel('Total Grating thickness (\mum)')
% ParamY3 = subplot(2,2,4);
% plot(Xparam,Yparam(:,3),'Marker','.')
% xlabel('Filling factor')
% ylabel('AR thickness (\mum)')


% %Total Grating thickness
% %-----------------------
% figure('FileName','Diamant')
% hold on
% nulldepth = subplot(1,3,1);
% title('AGPM Diamond/SiO2 : K-band')
% b=plot(Xparam,nulldpt,'Marker','o')
% xlabel('Total Grating thickness (\mum)')
% ylabel('Mean Null Depth')
% axis square
% set(gca,'YScale','log')
% ParamY1 = subplot(1,3,2);
% plot(Xparam,Yparam(:,1),'Marker','o')
% xlabel('Total Grating thickness (\mum)')
% ylabel('Period (\mum)')
% ParamY2 = subplot(1,3,3);
% plot(Xparam,Yparam(:,2),'Marker','o')
% xlabel('Total Grating thickness (\mum)')
% ylabel('AR thickness (\mum)')


%Sans balayage
%-------------
figure('FileName','Diamant')
hold on
%grid on
nulldepth = subplot(1,3,1);
plot(lb_t,nul_res_sp_b,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
axis square
set(nulldepth,'YScale','log')
xlabel('Wavelength \lambda (\mum)')
ylabel('Null Depth')
title('AGPM Diamond/SiO2 : K-band')
gtext({['Period \Lambda: ',num2str(Lb),' \mum'],['Total Grating thickness d: ',num2str(d),' \mum'],['AR thickness d_{AR}: ',num2str(d_AR),' \mum'],['Filling factor F: ',num2str(F*100),' %']})
transmittance = subplot(1,3,2);
tmp12=mean(tmp1+tmp2)/2;
plot(lb_t,tmp1,'-',lb_t,tmp2,'--',lb_t,tmp12,'k:','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('Transmittance')
legend('TE','TM')
phaseshift = subplot(1,3,3);
tmp3moy=mean(tmp3);
plot(lb_t,tmp3,lb_t,tmp3moy,'k:',lb_t,pi,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('TE-TM Phase Shift  \Delta\Phi_{TE-TM}  (rad)')