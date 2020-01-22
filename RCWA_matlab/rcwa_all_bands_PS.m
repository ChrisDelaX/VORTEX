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

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 dphi_sp_T 
global nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente

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

p=1;%13;%

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=15;%1;%30;%
pente=rad(0);
%pente=rad(4.3);
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

E_man=2.38;
E_man_min=1.499;%2.1;%
E_man_max=2.801;

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

%Pas du réseau (sublambda!)     --> limite Diamant bande N : 3.78
Lb=6;
Lb_min=0.1;Lb_max=7.6;

%Epaisseurs
d=12;
d_min=13;d_max=14;
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;

%double AGPM
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Facteurs de remplissage         
F=0.5;
F_min=0.56;F_max=0.65;

%Onde incidente 
%--------------

%Longueur d'onde
lb=10.5;
lb_min=10.5;lb_max=12.25;
nlb=41;%6;%21;%

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

disp('start')
datestr(now)
tic

% x0=[];
%x0=[Lb,F];

x0=[F,Lb,d];

%lb_vect=[1.15 1.5 2 3.45 4.35 8.5 11];
lb_vect=[3.45 4.35 8.5 11];
npts_lb=size(lb_vect,2);%8;%10;%20;%
npts_E=81;%11;%3;%

for i=1:npts_lb
    %x0
    i
    lb=lb_vect(i)
    for j=1:npts_E
        E_man=E_man_min+(E_man_max-E_man_min)*(j-1)/(npts_E-1);
        Lbzog = lb/E_man;
        x0(2)= 0.99*Lbzog;
        %x0(3)= 12;
        %dlim=0.8*Lbzog/(2*tan(pente))
        
        j
        %[x,fval,exitflag] = fminsearchbnd(@null,x0,[0.35 0.5*Lbzog 0],[0.8 Lbzog 7*Lbzog],optimset('MaxIter',15,'Display','iter','TolX',1e-2,'TolFun',1e-3));
        [x,fval,exitflag] = fminsearchbnd(@null_PS,x0,[0.01 0.5*Lbzog 0],[0.99 Lbzog inf],optimset('MaxIter',35,'Display','iter','TolX',5e-3,'TolFun',1e-5));
        %fval=null(x0);
        %x
        Xparam(j)=E_man;
        Yparam1(i,j)=x(1);
        Yparam2(i,j)=x(2);
        Yparam3(i,j)=x(3);
        
        esp_phi(i,j)=mean(tmp3);
        ecart_phi(i,j)=std(tmp3);
        esp_null(i,j)=mean(nul_res_sp_b);
        ecart_null(i,j)=std(nul_res_sp_b);
        
        nulldpt(i,j)=mean(nul_res_sp_b);
        rms_err(i,j)=fval;
        
    end
end




disp('end')
datestr(now)
toc

save all_spectrum_phase_shift_slope0_5



dessin_PS




