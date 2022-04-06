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
warning off MATLAB:illConditionedMatrix
warning off MATLAB:break_outside_of_loop
warning off MATLAB:declareGlobalBeforeUse


global WW sig prof prof_inv N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 limZOG p L pente d F E1 E2 

global N psi theta phi L
global E_arret_choice EI_choice EIII_choice E1_choice E2_choice E_AR_choice E_man


%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)  
N=8;%12;%20;%

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

p=13;%30;%3;%

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=15;%30;%
pente=rad(0);
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

% introduction manuelle !!! E_man = indice de réfraction
E_man=1.494;%1.44;
E_man_min=1.5;
E_man_max=2.5;

%Permittivités milieu extérieur incident   
EI_choice=1;

%Permittivités milieu extérieur émergent 
EIII_choice=8;

%Permittivités du réseau: EII
E1_choice=1;  %(la plus basse)
E2_choice=8;

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
Lb=1;%3.8;
Lb_min=4.7;Lb_max=4.75;

%Epaisseurs
d=3;%12.2;%8.87;%
d_min=10;d_max=14;
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage         
F=0.5;%0.327;%0.3547;
F_min=0.3;F_max=0.7;

%Onde incidente 
%--------------

%Longueur d'onde
lb=10;%0.55;%      35 9.5 738
lb_min=0.5;%0.3;%
lb_max=12;
nlb=41;%1;%

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

x0=[d,F];

disp('start')
datestr(now)
tic

npts_lb=3;%8;%10;%20;%
npts_E=3;%6;%5;%10;%

j=1;
for i=1:npts_lb
    i
    lb=lb_min+(lb_max-lb_min)*(i-1)/(npts_lb-1);
    Xparam(i)=lb;
    
    for j=1:npts_E
        j
        E_man=E_man_min+(E_man_max-E_man_min)*(j-1)/(npts_E-1);
        
%         x0=[d,F,Lb];
%         [x,fval,exitflag] = fminsearch(@null_ration,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-6,'TolFun',1e-6));
%         if x(3)>lb/E_man
%             Lb=lb/E_man;
%             x0=[d,F];
%             [x,fval,exitflag] = fminsearch(@null_ration_2,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-6,'TolFun',1e-6));
%         end
        
        Lb=lb/E_man;
        [x,fval,exitflag] = fminsearch(@null_ration_2,x0,optimset('MaxIter',15,'Display','iter','TolX',1e-6,'TolFun',1e-6));
        x
        
        Yparam(j,:)=E_man;
        nulldpt(i,j)=fval;
        Lbvar(i,j)=Lb;

    end
end

disp('end')
datestr(now)
toc

save ratio_all_spectrum

% Multi-param (lb,E_man)
%------------------

[vi,i]=min(nulldpt,[],1);
[min_nulldpt,j]=min(vi);

f=figure('FileName','Diamant');
set(f,'color',[1 1 1])
hold on
%grid on

hS = pcolor(Xparam,Yparam,log(nulldpt)');
shading interp
hC = colorbar;
%set(hC,'YScale','log')
set(hC,'FontSize',16,'FontWeight','bold')

[b]=plot(Xparam(i(j)),Yparam(j),'Marker','+');
set(b,'Linewidth',2,'color',[1 0 0]) 
xlabel('lb','FontSize',16,'FontWeight','bold')
ylabel('E','FontSize',16,'FontWeight','bold')
title('Null Depth','FontSize',16,'FontWeight','bold')
%set(gca,'YDir','reverse')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')



