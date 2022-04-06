%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Gestion de la m�moire
%---------------------
%---------------------

%Lib�ration
%----------

save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1p1 tmp1m1 tmp2p1 tmp2m1 tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof

%Param�tres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)
N=8;%12;%

%Allocation m�moire
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


%Lectures des entr�es
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
%13: profil trap�zoidal simple
%   -> L couches
%14: profil trap�zoidal avec couche antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2 couches
%15: profil trap�zoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trap�zoidal avec couche antireflet continue
%   -> L+1 couches
%17: profil trap�zoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trap�zoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2+1 couches
%19: profil trap�zoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trap�zoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trap�zoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trap�zoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adh�rent aux parois obliques
%   -> L+2+1 couches
%23: profil trap�zoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trap�zoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)

p=13;%1;%3;%

%Nombre de couches de discr�tisation du profil non rectangulaire en sus des milieux ext�rieurs
L=8;%12;%50;%

pente_deg=2.75;%0;%
pente=rad(pente_deg);


%Choix mat�riau --> Permittivit�s
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

E_man=1.5; % si manuel (8)

%Permittivit�s milieu ext�rieur incident
EI_choice=1;
%Permittivit�s milieu ext�rieur �mergent
EIII_choice=3;
%Permittivit�s du r�seau: EII
E1_choice=1;  %(la plus basse)
E2_choice=3;

%Permittivit�s de la couche antireflet
E_AR_choice=1;%14;%
E_AR1_choice=1;
E_AR2_choice=1;
%Permittivit�s de la couche d'arret
E_arret_choice=1;


%Param�tres du r�seau
%--------------------

%Pas du r�seau
Lb=1.42;%4.6;%
Lb_min=4.7;Lb_max=4.75;
%Epaisseurs
d=4.5;%13.69;%
d_min=3.5;%3.9;%
d_max=6.05;%7;%
d_AR=0.336;
d_AR_min=0.321;d_AR_max=0.351;
%Facteurs de remplissage
F=0.4;%0.2354;%0.19;%0.3547;%
F_min=0.35;%
F_max=0.51;%0.57;%

%tan(3.5/180*pi)*5.5*2+.52*1.42


%double AGPM
%-----------
%d=d/2;
%d_min=d_min/2;
%d_max=d_max/2;

%Onde incidente
%--------------

%Longueur d'onde
lb=3.5;%11;%
lb_min=3.4;%3.25;%
lb_max=4.1;%4.25;%
nlb=11;%21;%

%Angle incidence non conique
theta=rad(0);
theta_min=rad(30); theta_max=rad(50);
%Angle incidence conique
phi=rad(0);
phi_min=rad(0); phi_max=rad(0);
%Angle polarisation
psi=rad(45);
psi_min=rad(0); psi_max=rad(0);


% Filters
% -------

load filters.mat
x1 = filter_wide(:,1);
y1 = filter_wide(:,2);
p_wide = fit(x1,y1,'smoothingspline');


% Calcul Null / RMS_err  (optimistation)
%---------------------------------------

pente_min = 2.7;%2.5;%
pente_max = 3.2;%3;%
npts_pente = 3;%11;%5;%

npts_F=14;%55;%2;%5;%1440;%
npts_d=9;%33;%2;%5;%878;%

zoom='in';

x0=[];

disp('start')
tps0=now;
datestr(tps0)
tic

for k=1:npts_pente
    k
    pente_deg = pente_min+(pente_max-pente_min)*(k-1)/(npts_pente-1)
    pente=rad(pente_deg);
    
    for i=1:npts_F
        i;
        F=F_min+(F_max-F_min)*(i-1)/(npts_F-1);
        %F_equiv=F_min+(F_max-F_min)*(i-1)/(npts_F-1);
        %Lb=Lb_min+(Lb_max-Lb_min)*(i-1)/(npts_F-1);
        %d=d_min+(d_max-d_min)*(i-1)/(npts_F-1);
        
        Xparam(i)=F;
        
        if k==1 & i==2
            toc
            tpsfin = tps0 + npts_pente*npts_F*(now-tps0);
            disp('Fin pr�vue (date+heure) :')
            datestr(tpsfin)
        end
        
        for j=1:npts_d
            j;
            d=d_min+(d_max-d_min)*(j-1)/(npts_d-1);
            %F = F_equiv - d/Lb*tan(pente); %!!! et pas   F = F_equiv/(1+d/Lb*tan(pente));
            
            Yparam(j,:)=d;
            
            dmax_geom = (1-F)/2*Lb/tan(pente);
            if d < dmax_geom
                
                [fval1,fval2] = null(x0);
                %nuldpt(i,j) = fval1;
                %rms_err(i,j) = fval2;
                
                nul_tmp = nul_res_sp_b;
                nuldpt(i,j) = sum(nul_tmp.*p_wide(lb_t)')/sum(p_wide(lb_t));
                
            else
                nuldpt(i,j) = 1;
            end
            
        end
    end
    
    save(sprintf('multi_4L_zoom_%s_slope_%3.2f.mat',zoom,pente_deg))
    
    %multi_NULL_4L
end

disp('end')
datestr(now)
toc

multi_NULL_4L
%multi_RMS_L


