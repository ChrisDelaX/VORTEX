%Algorythme d'optimisation AGPM
%------------------------------
%------------------------------
%------------------------------


%Liberation
%----------
save backup.mat
clear all;close all;
warning off MATLAB:singularMatrix
warning off MATLAB:break_outside_of_loop
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 
global d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 
global E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice 
global E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente 
global fld1 n2lb limZOG prof
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost NARGghost 
global bandAR bandTAR nul_ghost nul_ghost_ARG nul_ideal
global tmp1p1 tmp1m1 tmp2p1 tmp2m1 optim



%Paramètres de calculs
%---------------------

%Troncature en X (N donc 2N+1 ordres au total)
N=6;%12;%
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
%14: profil trapézoidal avec couche antireflet dans les creux, sur les bosses 
%   et adhérent aux parois obliques
%   -> L+2 couches
%15: profil trapézoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trapézoidal avec couche antireflet continue
%   -> L+1 couches
%17: profil trapézoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, 
%    sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trapézoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches
%21: profil trapézoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, 
%    sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1 (sans couche d'arret)
%26: profil LETI 2 (avec couche d'arret)
p=13;%1;%3;%
%Nombre de couches de discrétisation du profil non rectangulaire en sus des milieux extérieurs
L=12;%16;%25;%50;%8;%


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
E_AR_choice=1;%14;
E_AR1_choice=1;
E_AR2_choice=1;
%Permittivités de la couche d'arret
E_arret_choice=1;


%Onde incidente
%--------------
%Angle incidence non conique
theta=deg2rad(0); theta_min=deg2rad(30); theta_max=deg2rad(50);
%Angle incidence conique
phi=deg2rad(0); phi_min=deg2rad(0); phi_max=deg2rad(0);
%Angle polarisation
psi=deg2rad(45); psi_min=deg2rad(0); psi_max=deg2rad(0);

% Select parameters
% ------------------
name = 'AGPM-METIS';
bands = string({'L','M','N1','N2'});%
params = string({'sw','lw', 'dep'});
nlb = 11;%8;%
nx = 32;
ny = 32;
nz = 1;
optim = 2;

%% Start simulation
% -----------------
tic
t0 = now;
t_start = datestr(t0)

% period limit for diamond at N: 3.78  L: 1.4277  K: 0.8391
% create_param: name, tag, units, val, delta, npts
first = 1;
for band = bands
for fixed_param = params
[band, fixed_param]
switch band
    case 'L'
%        lam = create_param('Wavelength', 'µm', 3.8, 0.6/3.8, nlb);
        lam = create_param('Wavelength', 'µm', 3.5, 1.2/3.8, nlb);
        per = create_param('Grating period', 'µm', 1.21, 0, 1);
        switch fixed_param
            case 'sw'
                lw = create_param('Line width', 'µm', 0.65, 0.1, nx);
                dep = create_param('Depth', 'µm', 5.53, 0.2, ny);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0, nz);
            case 'lw'
                lw = create_param('Line width', 'µm', 0.65, 0, nz);
                dep = create_param('Depth', 'µm', 5.53, 0.2, nx);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, ny);
            case 'dep'
                lw = create_param('Line width', 'µm', 0.65, 0.1, ny);
                dep = create_param('Depth', 'µm', 5.53, 0, nz);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, nx);
        end
    case 'M'
        lam = create_param('Wavelength', 'µm', 4.6, 1.4/4.6, nlb);
        per = create_param('Grating period', 'µm', 1.63, 0, 1);
        switch fixed_param
            case 'sw'
                lw = create_param('Line width', 'µm', 0.83, 0.1, nx);
                dep = create_param('Depth', 'µm', 6.55, 0.2, ny);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0, nz);
            case 'lw'
                lw = create_param('Line width', 'µm', 0.83, 0, nz);
                dep = create_param('Depth', 'µm', 6.55, 0.2, nx);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, ny);
            case 'dep'
                lw = create_param('Line width', 'µm', 0.83, 0.1, ny);
                dep = create_param('Depth', 'µm', 6.55, 0, nz);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, nx);
        end
    case 'N1'
        lam = create_param('Wavelength', 'µm', 9.25, 2.5/9.25, nlb);
        per = create_param('Grating period', 'µm', 3.36, 0.05, 1);
        switch fixed_param
            case 'sw'
                lw = create_param('Line width', 'µm', 1.64, 0.1, nx);
                dep = create_param('Depth', 'µm', 12.57, 0.2, ny);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0, nz);
            case 'lw'
                lw = create_param('Line width', 'µm', 1.64, 0, nz);
                dep = create_param('Depth', 'µm', 12.57, 0.2, nx);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, ny);
            case 'dep'
                lw = create_param('Line width', 'µm', 1.64, 0.1, ny);
                dep = create_param('Depth', 'µm', 12.57, 0, nz);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, nx);
        end
    case 'N2'
        lam = create_param('Wavelength', 'µm', 11.75, 3.5/11.75, nlb);
        per = create_param('Grating period', 'µm', 4.20, 0.05, 1);
        switch fixed_param
            case 'sw'
                lw = create_param('Line width', 'µm', 2.11, 0.1, nx);
                dep = create_param('Depth', 'µm', 16.57, 0.2, ny);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0, nz);
            case 'lw'
                lw = create_param('Line width', 'µm', 2.11, 0, nz);
                dep = create_param('Depth', 'µm', 16.57, 0.2, nx);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, ny);
            case 'dep'
                lw = create_param('Line width', 'µm', 2.11, 0.1, ny);
                dep = create_param('Depth', 'µm', 16.57, 0, nz);
                sw = create_param('Sidewall angle', 'deg', 2.45, 0.3, nx);
        end        
end
% AR grating, 3 values over 20% bandwidth
ARlam = create_param('AR wavelength', 'µm', lam.val, 0.2, 3);
ARtrans = [.993 .999 .995];

% assign x,y,z params
switch fixed_param
    case 'sw'
        xparam = lw;
        yparam = dep;
        zparam = sw;
    case 'lw'
        xparam = dep;
        yparam = sw;
        zparam = lw;
    case 'dep'
        xparam = sw;
        yparam = lw;
        zparam = dep;
end


%% Calculate Null Depth
% ---------------------
lb_t = lam.vals;
bandAR = ARlam.vals;
bandTAR = ARtrans;
Lb = per.val;
% material absorption
TH_diam_new
% start loops
for i=1:sw.npts %sidewall angle
    pente_deg = sw.vals(i);
    pente = deg2rad(pente_deg);    
    for j=1:lw.npts %line width
        j;
        F = lw.vals(j)/Lb;
        for k=1:dep.npts %depth
            d = dep.vals(k);
            [fval1,fval2] = null();
            if imag(fval1) == 0
                nulldepth = fval1;
            else
                nulldepth = 1;
            end
            switch fixed_param
                case 'sw'
                    nuldpt(j,k) = nulldepth;
                case 'lw'
                    nuldpt(k,i) = nulldepth;
                case 'dep'
                    nuldpt(i,j) = nulldepth;
            end
            if first==1
                first = 0;
                t_end = datestr(t0 + length(bands)*length(params)*...
                    toc/3600/24*sw.npts*lw.npts*dep.npts)
            end
        end
    end
end
Reject = 1./nuldpt;
dessin_multi
%fitswrite(nuldpt, sprintf('%s-%s_2D_%s=%3.2f.fits',name,band,fixed_param,zparam.val))
print('-dpng', sprintf('%s-%s_2D_%s=%3.2f.png',name,band,fixed_param,zparam.val), '-r300');
end
end
t_start
t_end
toc
datestr(now)
save('backup')