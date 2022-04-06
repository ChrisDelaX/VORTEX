%Gestion de la mémoire
%---------------------
%---------------------

    %Libération 
    %----------

    clear all;close all;

    %Lectures des entrées
    %--------------------
    %--------------------

    Inputs%_rough;%_old;
    
%     tt=round(rand(256));
%     
%     for i=1:1024,
%     t(i)=tt(floor(1+(i)*Lb/(1024*r)));
%     end;
%     for l=1:L_wo_R+2,
%         shift(l)=(cos(rand*pi)*r);
%     end;
    
    
    
    %Permittivités matériaux
    %-----------------------
    
    %1: Vide/Air
    %2: CdTe
    %3: Diamant
    %4: Germanium
    %5: Silicium
    %6: ZnSe
    %7: YF3
    
    Material;
                                   
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
    
%Scan1

        switch scan1
            case 1
                for ib=1:mesh1,
                Lb_t(ib)=Lbmin+ib*(Lbmax-Lbmin)/mesh1;
                Lb=Lb_t(ib);
                Scan2;
                saveHD;
                end;
            case 2
                for ib=1:mesh1,
                d_t(ib)=dmin+ib*(dmax-dmin)/mesh1;
                d=d_t(ib)
                Scan2;
                saveHD;
                end;
            case 3
                for ib=1:mesh1,
                Fa_t(ib)=Famin+ib*(Famax-Famin)/mesh1;
                Fa=Fa_t(ib);
                profile;%changement de profil => recalcul du profil
                Scan2;
                saveHD;
                end;
            case 4
                for ib=1:mesh1,
                Fb_t(ib)=Fbmin+ib*(Fbmax-Fbmin)/mesh1;
                Fb=Fb_t(ib);
                profile;%changement de profil => recalcul du profil
                Scan2;
                saveHD;
                end;
            case 5
                for ib=1:mesh1,
                theta_t(ib)=thetamin+ib*(thetamax-thetamin)/mesh1;
                theta=theta_t(ib);
                Scan2;
                saveHD;
                end;
            case 6
                for ib=1:mesh1,
                phi_t(ib)=phimin+ib*(phimax-phimin)/mesh1;
                phi=phi_t(ib);
                Scan2;
                saveHD;
                end;
            case 7
                for ib=1:mesh1,
                psi_t(ib)=psimin+ib*(psimax-psimin)/mesh1;
                psi=psi_t(ib);
                Scan2;
                saveHD;
                end;
            case 8
                for ib=1:mesh1,
                lb_t(ib)=lbmin+ib*(lbmax-lbmin)/mesh1;
                lb=lb_t(ib);
                Material;%changement de longueur d'onde => dispersion matériaux
                profile;%changement indice de profil => recalcul du profil
                Scan2;
%                 saveHD;
                end;
            case 9
                for ib=1:mesh1,
                d_AR_t(ib)=d_AR_min+ib*(d_AR_max-d_AR_min)/mesh1;
                d_AR=d_AR_t(ib);
                Scan2;
%                 saveHD;
                end;    
            otherwise
                ib=1;
                Scan2;
%                 saveHD;
        end;
            
%Résultats: Affichage
%--------------------

%Réflexion
%---------

% temp1(1:mesh5)=2*DER_s_fin(1,1,1,1,:,N+1);
% plot(lb_t,temp1.^2);hold on;
% temp2(1:mesh5)=2*DER_p_fin(1,1,1,1,:,N+1);
% plot(lb_t,temp2.^2);
% 
% temp3(1:mesh5)=unwrap(2*dphi_sp_R_fin(1,1,1,1,:));
% figure;plot(lb_t,temp3);
% q=(temp1./temp2).^2;
% 
% nul_phase_fin=(temp3-pi).^2;
% 
% nul_res_sp=((1+sqrt(q)).^2)./((1-sqrt(q)).^2+nul_phase_fin.*sqrt(q));
% 
% figure;plot(lb_t,1./nul_res_sp);

%Transmission
%------------

% temp1(1:mesh5)=2*DET_s_fin(1,1,1,1,:,N+1);
% %plot(lb_t,temp1);hold on;
% temp2(1:mesh5)=2*DET_p_fin(1,1,1,1,:,N+1);
% %plot(lb_t,temp2);
% 
% temp3(1:mesh5)=mod(dphi_sp_T_fin(1,1,1,1,:),2*pi);
% %plot(lb_t,temp3);
% mean(temp3)
% q=(temp1./temp2);
% 
% nul_phase_fin=(temp3-pi).^2;
% 
% nul_res_sp_1=((1+sqrt(q)).^2)./((1-sqrt(q)).^2+nul_phase_fin.*sqrt(q));
% 
% plot(lb_t,1./nul_res_sp_1);
