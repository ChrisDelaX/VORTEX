

clear all;
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
pente=rad(4);
VISIR_multi;
save VISIR4.mat

clear all;
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
pente=rad(4.5);
VISIR_multi;
save VISIR45.mat

clear all;
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
pente=rad(5);
VISIR_multi;
save VISIR5.mat

clear all;
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
pente=rad(5.5);
VISIR_multi;
save VISIR55.mat

clear all;
global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
pente=rad(6);
VISIR_multi;
save VISIR6.mat

break


% Fusion, mean, std
%------------------


clear all
load VISIR4.mat
nullmat4=nulldpt;
null_matrix=nullmat4;
VISIR_dessin
title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 4°','FontSize',16,'FontWeight','bold')

load VISIR45.mat
nullmat45=nulldpt;
null_matrix=null_matrix+nullmat45;
VISIR_dessin
title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 4.5°','FontSize',16,'FontWeight','bold')

load VISIR5.mat
nullmat5=nulldpt;
null_matrix=null_matrix+nullmat5;
VISIR_dessin
title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 5°','FontSize',16,'FontWeight','bold')

load VISIR55.mat
nullmat55=nulldpt;
null_matrix=null_matrix+nullmat55;
VISIR_dessin
title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 5.5°','FontSize',16,'FontWeight','bold')

load VISIR6.mat
nullmat6=nulldpt;
null_matrix=null_matrix+nullmat6;
VISIR_dessin
title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 6°','FontSize',16,'FontWeight','bold')

nulldpt=null_matrix./5;
VISIR_dessin
title('Mean Null Depth (3.5-4.1µm)  --  \Lambda = 1.42 µm / MEAN \alpha = 5°','FontSize',16,'FontWeight','bold')


% STD
%----


for i =1:size(nullmat4,1)
    for j=1:size(nullmat4,2)
        std_vect=[nullmat4(i,j) nullmat45(i,j) nullmat5(i,j) nullmat55(i,j) nullmat6(i,j)];
        std_matrix(i,j)=std(std_vect);
    end
end

nulldpt=std_matrix;
VISIR_dessin
title('Standard Deviation','FontSize',16,'FontWeight','bold')


