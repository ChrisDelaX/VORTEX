function[nul_res_sp_b]=null(x0)

global sig prof N lb lb_t nlb tmp1 tmp2 tmp3 tmp4 tmp5 dphi_sp_T 
global nul_res_sp_b null_res_sp prof
global ordre
global p d F Lb E1 E2

d=abs(x0(1));F=abs(x0(2));
Lb=abs(x0(3));


sig=ones(2*(2*N+1),nlb);

nlb=1;
for ib=1:nlb
    ib;
%    lb=lb_min+(lb_max-lb_min)*(ib-1)/(nlb-1);
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


