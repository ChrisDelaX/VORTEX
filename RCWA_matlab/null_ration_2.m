function[null_res_sp]=null_ration_2(x0)

global N lb lb_t nlb tmp1 tmp2 tmp3 tmp4 tmp5 dphi_sp_T 
global ordre
global p d L pente F Lb 

global N psi theta phi L
global E_arret_choice EI_choice EIII_choice E1_choice E2_choice E_AR_choice E_man

d=abs(x0(1));F=abs(x0(2));
%Lb=abs(x0(3));


for ib=1:nlb
    ib;
    
    lb_t(ib)=lb+(0.2*lb)*(ib-1)/(nlb-1);
    
    
    %Dispersion matériaux
    %--------------------

    [E1,E2,E_arret,E_AR,EI,EIII]=material(lb_t(ib));
        
    %Fabrication profile
    %-------------------

    [prof,d_t]=profile(p,L,d,F,pente,Lb,E1,E2,E_arret,E_AR);

    %Core RCWA principal
    %-------------------
    [DER_s,DER_p,DET_s,DET_p,dphi_sp_T,dphi_sp_R]=main(lb_t(ib),EI,EIII,Lb,prof,d_t);

%     if ib ==2
%         error('stop')
%     end
    
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


