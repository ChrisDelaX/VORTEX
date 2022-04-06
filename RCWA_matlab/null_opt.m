function[null_res_sp]=null(x0)

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 dphi_sp_T 
global nul_res_sp_b null_res_sp n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof
global ordre

%Lb=abs(x0(1));
F=abs(x0(1));d=abs(x0(2));d_AR=abs(x0(3));
%Lb=d/(1-F);
%F=abs(x0(2));
%Lb=abs(x0(3));

sig=ones(2*(2*N+1),nlb);


%for ib=1:nlb
for ib=1:1
    ib;
    %lb=lb_min+(lb_max-lb_min)*(ib-1)/(nlb-1);
    lb_t(ib)=lb;
    
    
    %Dispersion matériaux
    %--------------------
%     ib = 1
%     lb=10.5;
    material
%     n2lb(ib)=sqrt(E2)
%     limZOG(ib)=lb/n2lb(ib)
    
    
    %Optim AR
    %--------
%     EII=sqrt(EIII);
%     F_AR(ib)=(EII-EI)/(EIII-EI);
%     d_AR(ib)=lb/(4*sqrt(EII));
    

    %Fabrication profile
    %-------------------
    profile
    prof_inv=prof.^-1;
%     if ib == 1
%         figure; imagesc(prof)
%     elseif ib == nlb
%         figure; imagesc(prof)
%     end


    %Core RCWA principal
    %-------------------
    main
    sig(:,ib)=(Sigmatmp);
    % WW(:,:,ib)=W;
    %size(DET_s)
    %DET_s(N+1+ordre)
    ordre=0;
    temp1(ib)=2*DET_s(N+1+ordre);
    temp2(ib)=2*DET_p(N+1+ordre);
    temp4(ib)=2*DER_s(N+1+ordre);
    temp5(ib)=2*DER_p(N+1+ordre);
    
    % Simple AGPM
    temp3(ib)=mod(dphi_sp_T,2*pi);
    q(ib)=(temp1(ib)./temp2(ib));
    % Double AGPM
    % temp3(ib)=mod(2*dphi_sp_T,2*pi);
    % q(ib)=(temp1(ib)./temp2(ib))^2;
    
    nul_phase_fin(ib)=(temp3(ib)-pi).^2;
    
    %resultat du nulling avec phase et amplitude mismatch
    %----------------------------------------------------
    nul_res_sp_b(ib)=((1-sqrt(q(ib))).^2+nul_phase_fin(ib).*sqrt(q(ib)))./(1+sqrt(q(ib))).^2;

    %resultat du nulling avec phase mismatch only
    %--------------------------------------------
    % nul_res_sp_b(ib)=nul_phase_fin(ib)/4;
    
end

tmp1=(temp1);
tmp2=(temp2);
tmp3=(temp3);
tmp4=(temp4);
tmp5=(temp5);

null_res_sp=mean(nul_res_sp_b);


% % on ajoute l'aspect ratio comme critère d'optimisation
% if F < 0.5
%     F_aspRat = F;
% else
%     F_aspRat = 1-F;
% end
% aspRat = d/(F_aspRat*Lb);
% null_res_sp_AR = null_res_sp * aspRat^2;



% graphe optim AR
% ---------------
% figure
% plot(lb_t,F_AR)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Filling factor')
% figure
% plot(lb_t,d_AR)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Depth d (\mum)')





% null_res_sp=mean(nul_res_sp_b)+(mean(temp1(:).*0.01.*temp4(:).*0.99+temp2(:).*0.01.*temp5(:).*0.99))/2;%+(mean(temp4+temp5))/2;
% null_res_sp=(0.5*nul_res_sp_b(1)+nul_res_sp_b(8)+0.5*nul_res_sp_b(15)+0.5*nul_res_sp_b(21)+nul_res_sp_b(28)+0.5*nul_res_sp_b(32))/4;
% null_res_sp=mean(nul_res_sp_b(3))/2+mean(nul_res_sp_b(9))/2;
% 1/mean(nul_res_sp_b(4:37))
% 1/mean(nul_res_sp_b(64:95))
% null_res_sp=((mean(nul_res_sp_b(3:13)))); 
% null_res_sp=((mean(nul_res_sp_b(5:25))))%+mean(nul_res_sp_b(21:30))))/2;
% null_res_sp=((nul_res_sp_b(7))+(nul_res_sp_b(28)))/2;


% graphe Fabian
% -------------
% figure
% hold on
% plot(lb_t,fab1,lb_t,fab2);
% xlabel('Wavelength (microns)')
% ylabel('DET')
% legend('DET_s','DET_p',0)