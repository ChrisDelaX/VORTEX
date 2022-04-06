% Dessins rms error
% -----------------

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
lwz = 2;  % line width


pente_min = 2.5;
pente_max = 3;
npts_pente = 17;

npts_F=55;%2;%14;%1440;%
npts_d=33;%2;%9;%878;%


temp=0;
for i=1:npts_pente
    pente_deg = pente_min+(pente_max-pente_min)*(i-1)/(npts_pente-1);
    pente=rad(pente_deg);
    load(['multi_slope_',sprintf('%3.2f',pente_deg),'.mat'])
    temp = temp+rms_err;
    
    %steps
    %-----
    fmulti = rms_err;
    dessin_multi
    caxis([-2,0.5]);  %large
    %caxis([-1.7,0]);  %zoom
    t=title(['                    RMS Phase Shift Error: log_{10}( \epsilon_{rms}) for \alpha = ',sprintf('%3.2f',pente_deg),'°']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
    set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
    saveas(gcf,['RMS_multi_slope_',sprintf('%3.2f',pente_deg),'.fig']);
    print('-dpng', ['RMS_multi_slope_',sprintf('%3.2f',pente_deg),'.png'], '-r300'); 
end
mu_nu = temp./npts_pente;


% mean
% ----
fmulti = mu_nu;
dessin_multi
caxis([-2,0.5]);  %large
%caxis([-1.7,0]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
t=title(['Expected Value (mean) of \epsilon_{rms} (log_{10})']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
saveas(gcf,['RMS_multi_MEAN.fig']);
print('-dpng', ['RMS_multi_MEAN.png'], '-r300');



temp=0;
for i=1:npts_pente
    pente_deg = pente_min+(pente_max-pente_min)*(i-1)/(npts_pente-1);
    pente=rad(pente_deg);
    load(['multi_slope_',sprintf('%3.2f',pente_deg),'.mat'])
    temp = temp+(rms_err-mu_nu).^2;
end
s_nu = sqrt(temp./(npts_pente-1));

% std
% ---
fmulti = s_nu;
dessin_multi
caxis([-4,-1]);  %large
%caxis([-4,-1]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
t=title(['Standard Deviation of \epsilon_{rms} (log_{10})']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
saveas(gcf,['RMS_multi_STD.fig']);
print('-dpng', ['RMS_multi_STD.png'], '-r300');

