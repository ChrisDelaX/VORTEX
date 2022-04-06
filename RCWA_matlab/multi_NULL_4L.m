% Dessins mean null depth
% -----------------------

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
lwz = 1;%2;  % line width


pente_min = 2.7;%2.5;%
pente_max = 3.2;%3;%
npts_pente = 11;%3;%5;%

npts_F=55;%14;%2;%5;%1440;%
npts_d=33;%9;%2;%5;%878;%

zoom='in';


pente_minz = pente_min;
pente_maxz = pente_max;
npts_pentez = npts_pente;
temp=0;
for i=1:npts_pentez
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(i-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    load(sprintf('multi_4L_zoom_%s_slope_%3.2f.mat',zoom,pente_deg))
    temp = temp+nuldpt;
    
    %steps
    %-----
    fmulti = nuldpt;
    dessin_multi
    %dessin_multi_normal
    caxis([-3.6,-0.4]);  %large
    set(hC,'ytick',[-3.5 -3 -2.5 -2 -1.5 -1 -.5])
    %caxis([-1.7,0]);  %zoom
    t=title(['Mean Null Depth: log_{10}(N) for \alpha = ',sprintf('%3.2f',pente_deg),'°']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
    set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
    %set(t,'FontSize',12,'HorizontalAlignment','center')
    grid on
%    print('-dpng',sprintf('NULL_multi_4L_zoom_%s_slope_%3.2f.png',zoom,pente_deg), '-r300')
%    print('-depsc2',sprintf('NULL_multi_4L_zoom_%s_slope_%3.2f.eps',zoom,pente_deg), '-r300')
end
mu_nu = temp./npts_pentez;


% mean
% ----
fmulti = mu_nu;
dessin_multi
%dessin_multi_normal
caxis([-3.1,-0.4]);  %large
set(hC,'ytick',[-3 -2.5 -2 -1.5 -1 -.5])
%caxis([-1.7,0]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
%t=title(['Expected Value (mean) of Mean Null Depth (log_{10})']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
t=title(['Mean Null Depth (log_{10}) for \alpha = ',sprintf('%3.2f',pente_minz),'-',sprintf('%3.2f',pente_maxz),'°']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
%set(t,'FontSize',12,'HorizontalAlignment','center')
%print('-dpng',sprintf('NULL_multi_4L_zoom_%s_MEAN.png',zoom), '-r300')
print('-depsc2',sprintf('NULL_multi_4L_zoom_%s_MEAN.eps',zoom), '-r300')


temp=0;
for i=1:npts_pentez
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(i-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    load(sprintf('multi_4L_zoom_%s_slope_%3.2f.mat',zoom,pente_deg))
    temp = temp+(nuldpt-mu_nu).^2;
end
s_nu = sqrt(temp./(npts_pentez-1));

% std
% ---
fmulti = s_nu;
dessin_multi
%dessin_multi_normal
caxis([-3.8,-.7]);  %large
set(hC,'ytick',[-3.5 -3 -2.5 -2 -1.5 -1])
%caxis([-4,-1]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
t=title(['Standard Deviation (log_{10}) for \alpha = ',sprintf('%3.2f',pente_minz),'-',sprintf('%3.2f',pente_maxz),'°']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' µm] and \Lambda = ',sprintf('%3.2f',Lb),' µm']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
%set(t,'FontSize',12,'HorizontalAlignment','center')
%print('-dpng',sprintf('NULL_multi_4L_zoom_%s_STD.png',zoom), '-r300')
print('-depsc2',sprintf('NULL_multi_4L_zoom_%s_STD.eps',zoom), '-r300')

