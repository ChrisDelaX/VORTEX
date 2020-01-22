
function[null_res_sp,RMSErr]=null(x0)
%function[NTOT,null_res_sp]=null(x0)
%function[RMSErr,null_res_sp]=null(x0)
%function[Tin,null_res_sp]=null(x0)
%function[meritfun,null_res_sp]=null(x0)

global WW sig prof N p L lb lb_t nlb pente Lb d d_AR d_AR1 d_AR2 
global d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 
global E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice 
global E_arret E_arret_choice theta phi psi tmp1p1 tmp1m1 tmp2p1 tmp2m1 
global tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 dphi_sp_T dphi_sp_R
%global nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew fld1 n2lb 
global limZOG prof ghost lbtmp retardghost dn_sp_T n_s_T n_p_T LFt
global Tin Rin T0 R0 absor TARG RARG Ttot Rtot TARGtot RARGtot Nghost 
global NARGghost bandAR bandTAR nul_ghost nul_ghost_ARG nul_ideal 
global optim


switch optim
    case 0
        d=abs(x0(1));
    case 1
        F=abs(x0(1));d=abs(x0(2));
end
%[F*Lb d]


sig=ones(2*(2*N+1),nlb);

%lb_t = lbtmp ./1000;

SSErr = 0;
for ib=1:size(lb_t,2)
    
    ib;

    lb=lb_t(ib);
    
    %Dispersion mat�riaux
    %--------------------
    material
    n2lb(ib)=sqrt(E2);
    limZOG(ib)=lb/n2lb(ib);
    

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
    temp1(ib)=2*DET_s(N+1);
    temp2(ib)=2*DET_p(N+1);
    temp4(ib)=2*DER_s(N+1);
    temp5(ib)=2*DER_p(N+1);
    
    tmp1p1(ib)=2*DET_s(N+2);
    tmp1m1(ib)=2*DET_p(N+2);
    tmp2p1(ib)=2*DET_s(N);
    tmp2m1(ib)=2*DET_p(N);
    
    % Simple AGPM
    % -----------
    temp3(ib)=mod(dphi_sp_T,2*pi);
    retard(ib) = temp3(ib)/(2*pi);
    temp6(ib)=mod(dphi_sp_R,2*pi);
    retardghost(ib)= mod(temp3(ib)+temp6(ib),2*pi)/(2*pi);
 
%     % Double AGPM
%     % -----------
%     temp3(ib)=mod(2*dphi_sp_T,2*pi);
%     q(ib)=(temp1(ib)./temp2(ib))^2;
   
    
%     % dn trans
%     % --------
%     n_s_T(ib) = mod(phi_s_T,2*pi)*lb/(2*pi*d);
%     n_p_T(ib) = mod(phi_p_T,2*pi)*lb/(2*pi*d);
%     dn_sp_T(ib) = abs(n_s_T(ib)-n_p_T(ib));
    
    
    %resultat du nulling avec phase et amplitude mismatch
    %----------------------------------------------------
    nul_phase_fin(ib)=(temp3(ib)-pi).^2;
    q(ib)=(temp1(ib)./temp2(ib));
    nul_res_sp_b(ib)=((1-sqrt(q(ib))).^2+nul_phase_fin(ib).*sqrt(q(ib)))./(1+sqrt(q(ib))).^2;

%     %resultat du nulling avec phase mismatch only
%     %--------------------------------------------
%     nul_res_sp_b(ib)=nul_phase_fin(ib)/4;
    
    SSErr = SSErr + (temp3(ib)-pi)^2;

% null_res_sp=mean(nul_res_sp_b)+(mean(temp1(:).*0.01.*temp4(:).*0.99+ ...
%   temp2(:).*0.01.*temp5(:).*0.99))/2;%+(mean(temp4+temp5))/2;

end

tmp1=(temp1);
tmp2=(temp2);
tmp3=(temp3);
tmp4=(temp4);
tmp5=(temp5);
tmp6=(temp6);

% % FOR METIS
% null_res_sp = mean(nul_res_sp_b(nlb+1:end));
% phase = (mod(tmp3(1:nlb)+pi,2*pi)-pi)./pi;
% phase(isnan(phase)) = 1;
% meritfun = null_res_sp * std(phase);
% return

null_res_sp = mean(nul_res_sp_b);

MSErr = SSErr/(nlb-1);
RMSErr = sqrt(MSErr);


%%%

%in:AGPM
%Tin =pchip(lb_t,(tmp1+tmp2)./2,lb_t);
%Rin =pchip(lb_t,(tmp4+tmp5)./2,lb_t);
Tin =(tmp1+tmp2)./2;
Rin =(tmp4+tmp5)./2;

%out:

% WITHOUT ARG
% -----------
TNUM = T0.*Tin.*absor;
TDEN = 1-R0.*Rin.*absor.^2;
Ttot = TNUM ./ TDEN;
RNUM = R0.*Tin.^2.*absor.^2;
RDEN = 1-R0.*Rin.*absor.^2;
Rtot = Rin + RNUM ./ RDEN;
TnoARG=mean(Ttot);
RnoARG=mean(Rtot);

% WITH ARG
% --------
TARG=polyval(polyfit(bandAR,bandTAR,2),lb_t);
TARG(TARG<T0)=T0(TARG<T0);
RARG=1-TARG;

TNUM = TARG.*Tin.*absor;
TDEN = 1-RARG.*Rin.*absor.^2;
TARGtot = TNUM ./ TDEN;
RNUM = RARG.*Tin.^2.*absor.^2;
RDEN = 1-RARG.*Rin.*absor.^2;
RARGtot = Rin + RNUM ./ RDEN;
TwithARG=mean(TARGtot);
RwithARG=mean(RARGtot);

% GHOST
% -----
Tghost = Ttot-T0.*Tin.*absor;
Nghost = Tghost./Ttot;
TARGghost = TARGtot-TARG.*Tin.*absor;
NARGghost = TARGghost./TARGtot;
NTOT=mean(nul_res_sp_b+NARGghost);

% NULLS
% -----
nul_ghost = nul_res_sp_b+Nghost;
nul_ghost_ARG = nul_res_sp_b+NARGghost;
nul_ideal = nul_res_sp_b;
