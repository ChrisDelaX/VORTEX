T=298;    

% E_arret_choice=1;

switch E_arret_choice
    case 1
        E_arret=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        E_arret=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        lb=3.3%3.426%1.95;%1.48;%.6328;%1.15;%3.4;%
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        E_arret=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
        sqrt(E_arret)
        Lb=lb/sqrt(E_arret)
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        E_arret=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        E_arret=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        E_arret=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        E_arret=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        E_arret=E_man^2;  
    case 9
        asga_data;
        E_arret=(spline(asga_l,asga_n,lb)-sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        E_arret=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        E_arret=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        E_arret=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        E_arret=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        %lb=.882;%.980;%.430;%.550;%
        infrasil_data;
        E_arret=(spline(infra_l,infra_n,lb))^2;
        %sqrt(E_arret)
        %Lb=lb/sqrt(E_arret)
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        E_arret=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
    case 16
        si3n4;
        E_arret=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        %             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        %             E_arret_tmp=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
        A_zns=5.608e-5*T + 2.282; B_zns = -8.671e-6*T - 1.563e-2;C_zns = 5.549e-7*T + 2.067e-3; D_zns = 2.597e-8*T - 1.714e-4;E_zns = -9.798e-10*T + 2.884e-6;
        E_arret=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
        
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        E_arret=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        E_arret=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        E_arret=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

switch EI_choice
    case 1
        EI=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        EI=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        EI=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        EI=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        EI=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        EI=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        EI=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        EI=E_man^2;  
    case 9
        asga_data;
        EI=(spline(asga_l,asga_n,lb)-sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        EI=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        EI=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        EI=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        EI=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        infrasil_data;
        EI=(spline(infra_l,infra_n,lb))^2;
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        EI=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
    case 16
        si3n4;
        EI=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        %             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        A_zns=5.608e-5*T + 2.282; B_zns = -8.671e-6*T - 1.563e-2;C_zns = 5.549e-7*T + 2.067e-3; D_zns = 2.597e-8*T - 1.714e-4;E_zns = -9.798e-10*T + 2.884e-6;
        EI=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        EI=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        EI=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        EI=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

switch EIII_choice
    case 1
        EIII=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        EIII=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        EIII=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        EIII=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        EIII=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        EIII=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        EIII=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        EIII=E_man^2;  
    case 9
        asga_data;
        EIII=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        EIII=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        EIII=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        EIII=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        EIII=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        infrasil_data;
        EIII=(spline(infra_l,infra_n,lb))^2;
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        EIII=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
    case 16
        si3n4;
        EIII=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        EIII=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        EIII=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        EIII=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);     
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        EIII=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

switch E1_choice
    case 1
        E1=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        E1=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        E1=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        E1=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        E1=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        E1=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        E1=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        E1=E_man^2;
    case 9
        asga_data;
        E1=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        E1=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        E1=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        E1=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        E1=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        infrasil_data;
        E1=(spline(infra_l,infra_n,lb))^2;
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        E1=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
    case 16
        si3n4;
        E1=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        %             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        A_zns=5.608e-5*T + 2.282; B_zns = -8.671e-6*T - 1.563e-2;C_zns = 5.549e-7*T + 2.067e-3; D_zns = 2.597e-8*T - 1.714e-4;E_zns = -9.798e-10*T + 2.884e-6;
        
        E1=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        E1=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        E1=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        E1=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

switch E2_choice
    case 1
        E2=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        E2=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        E2=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        E2=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        E2=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        E2=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        E2=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        E2=E_man^2;  
    case 9
        asga_data;
        E2=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;
        
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        E2=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        E2=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        E2=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        E2=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        infrasil_data;
        E2=(spline(infra_l,infra_n,lb))^2;
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        E2=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);  
    case 16
        si3n4;
        E2=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        %             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        A_zns=5.608e-5*T + 2.282; B_zns = -8.671e-6*T - 1.563e-2;C_zns = 5.549e-7*T + 2.067e-3; D_zns = 2.597e-8*T - 1.714e-4;E_zns = -9.798e-10*T + 2.884e-6;
        
        E2=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        E2=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        E2=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        E2=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

switch E_AR_choice
    case 1
        E_AR=1;
    case 2
        T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
        E_AR=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
    case 3
        A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
        E_AR=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
    case 4
        T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
        E_AR=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
    case 5
        T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
        E_AR=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
    case 6
        T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
        E_AR=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
    case 7
        YF3_data;
        E_AR=(spline(YF3_l,YF3_n,lb)+sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
    case 8
        E_AR=E_man^2;
    case 9
        asga_data;
        E_AR=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;            
    case 10
        K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
        K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
        E_AR=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
    case 11
        GASIR_2;
        E_AR=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 12
        GASIR_1;
        E_AR=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
    case 13
        inp;
        E_AR=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
    case 14
        infrasil_data;
        E_AR=(spline(infra_l,infra_n,lb))^2;
    case 15
        A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
        E_AR=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
    case 16
        si3n4;
        E_AR=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
    case 17
        % A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
        A_zns=5.608e-5*T + 2.282; B_zns = -8.671e-6*T - 1.563e-2;C_zns = 5.549e-7*T + 2.067e-3; D_zns = 2.597e-8*T - 1.714e-4;E_zns = -9.798e-10*T + 2.884e-6;
        
        E_AR=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
    case 18
        K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
        K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
        E_AR=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
    case 19
        A_znseb=1.509e-4*T + 2.407; B_znseb = -1.801e-5*T - 2.564e-4;C_znseb = 1.300e-6*T - 1.308e-5 ; D_znseb = -3.878e-8*T - 1.480e-5;
        E_AR=sellmeier_znseb(lb,A_znseb,B_znseb,C_znseb,D_znseb);
    case 20
        A1_LAH83=1.82331579; A2_LAH83=5.48625885e-1 ; A3_LAH83= 1.63182855;B1_LAH83=9.11468875e-3;B2_LAH83=3.2841967e-2;B3_LAH83=1.23611174e+2;
        E_AR=sellmeier_lah83(lb,A1_LAH83,B1_LAH83,A2_LAH83,B2_LAH83,A3_LAH83,B3_LAH83);
end

%     switch E_AR1_choice
%         case 1
%             E_AR1=1;
%         case 2
%             T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
%             E_AR1=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
%         case 3
%             A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
%             E_AR1=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
%         case 4
%             T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
%             E_AR1=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
%         case 5
%             T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
%             E_AR1=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
%         case 6
%             T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
%             E_AR1=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
%         case 7
%             YF3_data;
%             E_AR1=(spline(YF3_l,YF3_n,lb)-sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
%         case 8
%             E_AR1=E_man^2;
%         case 9
%             asga_data;
%             E_AR1=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;            
%         case 10
%             K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
%             K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
%             E_AR1=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
%         case 11
%             GASIR_2;
%             E_AR1=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
%         case 12
%             GASIR_1;
%             E_AR1=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
%         case 13
%             inp;
%             E_AR1=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
%        case 14
%             infrasil_data;
%             E_AR1=(spline(infra_l,infra_n,lb))^2;
%        case 15
%             A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
%             E_AR1=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
%        case 16
%             si3n4;
%             E_AR1=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
%        case 17
%             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
%             E_AR1=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
%        case 18
%             K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
%             K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
%             E_AR1=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
%     end
%     
%     switch E_AR2_choice
%         case 1
%             E_AR2=1;
%         case 2
%             T_cdte=T; A_cdte = -2.973e-4*T_cdte + 3.8466; B_cdte = 8.057e-4*T_cdte + 3.2215; C_cdte = -1.10e-4*T_cdte + 0.1866; D_cdte = -2.160e-2*T_cdte + 12.718; E_cdte = -3.160e1*T_cdte + 18753;
%             E_AR2=sellmeier_cdte_ge(lb,A_cdte,B_cdte,C_cdte,D_cdte,E_cdte);
%         case 3
%             A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ; 
%             E_AR2=sellmeier_diamant(lb,A_d,B_d,C_d,D_d,E_d);
%         case 4
%             T_Ge=T; A_Ge=-6.040e-3*T_Ge + 11.05128 ; B_Ge = 9.295e-3*T_Ge + 4.00536 ; C_Ge = -5.392e-4*T_Ge + 0.599034 ; D_Ge = 4.151e-4*T_Ge + 0.09145; E_Ge = 1.51408*T_Ge + 3426.5 ;
%             E_AR2=sellmeier_cdte_ge(lb,A_Ge,B_Ge,C_Ge,D_Ge,E_Ge);;
%         case 5
%             T_si=T; A_si=1.600e-4*T_si+3.431; B_si=-2.643e-2; C_si=4.324e-3; D_si=-3.194e-4; E_si=8.835e-6;
%             E_AR2=sellmeier_Si(lb,A_si,B_si,C_si,D_si,E_si);
%         case 6
%             T_ZnSe=T; A_ZnSe = 1; B_ZnSe = 4.46395; C_ZnSe = 0.20108^2; D_ZnSe = 0.46132; E_ZnSe = 0.39211^2; F_ZnSe = 2.88289 ; G_ZnSe = 47.04759^2 ;
%             E_AR2=sellmeier_ZnSe(lb,A_ZnSe,B_ZnSe,C_ZnSe,D_ZnSe,E_ZnSe,F_ZnSe,G_ZnSe);
%         case 7
%             YF3_data;
%             E_AR2=(spline(YF3_l,YF3_n,lb)-sqrt(-1)*(spline(YF3_l,YF3_k,lb)))^2;
%         case 8
%             E_AR2=E_man^2;
%         case 9
%             asga_data;
%             E_AR2=(spline(asga_l,asga_n,lb)+sqrt(-1)*(spline(asga_l,asga_k,lb)))^2;            
%         case 10
%             K1_n_laf32=1.91731106;L1_n_laf32=9.78600358e-3;K2_n_laf32=2.24019825e-1;L2_n_laf32=3.86141473e-2;
%             K3_n_laf32=1.22087075;L3_n_laf32=8.46431265e1;
%             E_AR2=sellmeier_n_laf32(lb,K1_n_laf32,L1_n_laf32,K2_n_laf32,L2_n_laf32,K3_n_laf32,L3_n_laf32);
%         case 11
%             GASIR_2;
%             E_AR2=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
%         case 12
%             GASIR_1;
%             E_AR2=(spline(GAS2_l,GAS2_n,lb)+sqrt(-1)*(spline(GAS2_l,GAS2_k,lb)))^2;
%         case 13
%             inp;
%             E_AR2=(spline(l_inp,n_inp,lb)+sqrt(-1)*(spline(l_inp,k_inp,lb)))^2;
%        case 14
%             infrasil_data;
%             E_AR2=(spline(infra_l,infra_n,lb))^2;
%        case 15
%             A1_KRS5=3.75415165;B1_KRS5=2.08331008e-1;A2_KRS5=9.09002797e-1;B2_KRS5=3.77160136e-1;A3_KRS5=1.25434711e+1;B3_KRS5=1.65646424e+2;
%             E_AR2=sellmeier_krs5(lb,A1_KRS5,B1_KRS5,A2_KRS5,B2_KRS5,A3_KRS5,B3_KRS5);
%        case 16
%             si3n4;
%             E_AR2=(spline(si3n4_l,si3n4_n,lb)-sqrt(-1)*(spline(si3n4_l,si3n4_k,lb)))^2;
%        case 17
%             A_zns=2.29819;B_zns=-1.798e-2;C_zns=2.19e-3;D_zns=-1.614e-4;E_zns=2.538e-6;
%             E_AR2=sellmeier_Si(lb,A_zns,B_zns,C_zns,D_zns,E_zns);
%        case 18
%             K1_n_lasf44=1.78897105;L1_n_lasf44=8.725062277e-3;K2_n_lasf44=3.8675867e-1;L2_n_lasf44=3.08085023e-2;
%             K3_n_lasf44=1.30506243;L3_n_lasf44=9.27743824e1;
%             E_AR2=sellmeier_n_laf32(lb,K1_n_lasf44,L1_n_lasf44,K2_n_lasf44,L2_n_lasf44,K3_n_lasf44,L3_n_lasf44);
%     end
