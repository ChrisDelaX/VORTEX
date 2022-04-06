%Initialisation de l'algo d'optimisation
%---------------------------------------

%Gestion de la mémoire
%---------------------
%---------------------

%Libération 
%----------

clear all
close all
warning off MATLAB:singularMatrix

global WW sig prof N Lb_tmp d_tmp d_AR_tmp d_AR1_tmp d_AR2_tmp d_arret_tmp F_tmp E_man EI_choice EI EIII_choice EIII E1_choice E2_choice E_AR_choice E_AR1_choice E_AR2_choice E_arret_choice p L_tmp lb_t lb lbmin lbmax nlb theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 

%Lectures des entrées
%--------------------
%--------------------

%Inputs
%------

%Paramètres de calculs
%---------------------

%Troncature
%En X (N donc 2N+1 ordres au total)  

N=0;

%Paramètres du réseau
%--------------------    


%Pas du réseau
Lb_tmp=1.5;

%Epaisseur totale
d_tmp=1.2000;
dmin=12;dmax=22;
d_AR_tmp=0.35;
d_AR1_tmp=0.1185;
d_AR2_tmp=0.0715;
d_arret_tmp=0.1;

%Facteurs de remplissage moyen         
F_tmp=0.55;Famin=0.3;Famax=0.7;
%F_tmp=1-0.112/Lb_tmp;

%Choix du profil

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

p=1;

%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L_tmp=1;


%Permittivités
%-------------

%Choix matériau

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

E_man=1.46;
%LaF2: E_man=1.57;
%BaF2: E_man(l=1.97 microns)=1.4647
%BaF2: E_man(l=2.3253 microns)=1.4635
%BaF2: E_man(l=2.6738 microns)=1.4623

%Permittivités milieu extérieur incident   
EI_choice=1;

%Permittivités milieu extérieur émergent 
EIII_choice=14;

%Permittivités du réseau
%-----------------------

%Permittivité milieu 1 (la plus basse)
E1_choice=1;

%Permittivité milieu 2
E2_choice=14;

%Permittivités de la couche antireflet
%-------------------------------------

E_AR_choice=14;

E_AR1_choice=7; E_AR2_choice=17;

%Permittivités de la couche d'arret
%--------------------------------

E_arret_choice=16;


%Onde incidente 
%--------------

%Longueur d'onde
lb=4;
lbmin=1.99;lbmax=2.3;
nlb=16; %nombre de bandes de lb

%Angle incidence non conique
theta= 0*pi/180	;
thetamin=30*pi/180;thetamax=50*pi/180;

%Angle incidence conique
phi=0*pi/180;
phimin=0;phimax=0;

%Angle polarisation
psi=pi/4;
psimin=0;psimax=0;


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


%Optimisation
%------------
%     mfz=F_tmp*Lb_tmp;
% maxfz=Lb_tmp*F_tmp;
% minfz=Lb_tmp*(1-F_tmp);
% x0=[maxfz,minfz,d_tmp,d_AR_tmp]
%     x0=[Lb_tmp,d_tmp,F_tmp,d_AR_tmp];
%     x0=[Lb_tmp,F_tmp,d_tmp,d_AR_tmp]
x0=[Lb_tmp,F_tmp];

for i=1:20,
    %       %  %i=1;
    d(i)=dmin+(dmax-dmin)*i/19;
    
    d_tmp=d(i);
    
    %         Lb_t(i)=0.33+(i-1)*0.01;
    %         Lb_tmp=Lb_t(i);
    %         F_t(i)=1-0.08/Lb_t(i);
    %         F_tmp=F_t(i);
    
    
    %         Fa_t(i)=Famin+(i-1)*(Famax-Famin)/9;F_tmp=Fa_t(i);
    %x0=x0_t(i,:);
    %Fa=Fa_t(i);
    [x,fval,exitflag] = fminsearch(@null_new,x0,optimset('MaxIter',50,'Display','iter','TolX',5e-3,'TolFun',1e-5))
    %          parspace=[0.3 0.45; 0.1 0.9; 1 5; 0.05 0.95];
    %          opts=[1e-6 1 1];
    %         options=[20,-1,0.12,10,20000,2000,1e-6,50,0];
    %         x=ga(parspace,'null_new');
    %         [x,fval,exitflag] = fminsearch(@null_new,x0,optimset('MaxIter',200,'Display','iter','TolX',1e-4,'TolFun',1e-10))
    %   [x,endPop,bPop,traceInfo]= ga(parspace,'null_new',[],[],opts)      
    %[x,fval,exitflag] = fminimax(@null_tirg,x0,[],[],[],[],[],[],[],optimset('MaxIter',200,'Display','iter'))
    %[x,fval,exitflag] = fgoalattain(@null_tirg,x0,1e-9,2,[],[],[],[],[],[],[],optimset('MaxIter',200,'Display','iter'))
    
    param(i,:)=x;
    nulldpt(i)=fval;
    %         tm1(i,:)=tmp1;
    %         tm2(i,:)=tmp2;
    %         tm3(i,:)=tmp3;
    %         tm4(i,:)=tmp4;
    %         tm5(i,:)=tmp5;
    %       %  
end;
%         
%         end;



