% Dessins mean null depth
% -----------------------

global WW sig prof N p L lb lb_min lb_max lb_t nlb Lb d d_AR d_AR1 d_AR2 d_arret d_t F E_man EI EI_choice EIII EIII_choice E1 E1_choice E2 E2_choice E_AR E_AR_choice E_AR1 E_AR1_choice E_AR2 E_AR2_choice E_arret E_arret_choice theta phi psi tmp1 tmp2 tmp3 tmp4 tmp5 
global dphi_sp_T nul_res_sp_b null_res_sp retard n1 n2 n3 Fnew dnew pente fld1 n2lb limZOG prof

fnz = 'Arial'; % fontname
fsz = 28; % fontsize
fwz = 'normal';%'Bold'; % fontweight
lwz = 2;  % line width


pente_min = 2.7;%2.5;%
pente_max = 3.2;%3;%
npts_pente = 3;%11;%5;%

npts_F=14;%55;%2;%9;%1440;%
npts_d=9;%33;%2;%5;%878;%

pente_minz = pente_min;
pente_maxz = pente_max;
npts_pentez = npts_pente;
temp=0;
for i=1:npts_pentez
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(i-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    load(['multi_slope_L_',sprintf('%3.2f',pente_deg),'.mat'])
    temp = temp+nuldpt;
    
    %steps
    %-----
    fmulti = nuldpt;
    dessin_multi
    %set(gca,'xlim',[.3 .57])
    set(gca,'xtick',[.3 .35 .4 .45 .5])
    %set(gca,'ytick',[4 4.25 4.5 4.75 5 5.25])
    %set(gca,'ylim',[4 5.25])
    caxis([-3.8,-0.2]);  %large
    set(hC,'ytick',[-3.5 -3 -2.5 -2 -1.5 -1 -.5])
    %caxis([-1.7,0]);  %zoom
    t=title(['                    Mean Null Depth: log_{10}(N) for \alpha = ',sprintf('%3.2f',pente_deg),'�']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' �m] and \Lambda = ',sprintf('%3.2f',Lb),' �m']);
    set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
    %saveas(gcf,['NULL_multi_slope_L_',sprintf('%3.2f',pente_deg),'.fig']);
    print('-dpng', ['NULL_multi_slope_L_',sprintf('%3.2f',pente_deg),'.png'], '-r300'); 
end
mu_nu = temp./npts_pentez;


% mean
% ----
fmulti = mu_nu;
dessin_multi
%set(gca,'xlim',[.3 .57])
set(gca,'xtick',[.3 .35 .4 .45 .5])
%set(gca,'ytick',[4 4.25 4.5 4.75 5 5.25])
%set(gca,'ylim',[4 5.25])
caxis([-3.3,-0.2]);  %large
set(hC,'ytick',[-3 -2.5 -2 -1.5 -1 -.5])
%caxis([-1.7,0]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
%t=title(['Expected Value (mean) of Mean Null Depth (log_{10})']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' �m] and \Lambda = ',sprintf('%3.2f',Lb),' �m']);
t=title(['     Mean Null Depth (log_{10}) for \alpha = ',sprintf('%3.2f',pente_minz),'-',sprintf('%3.2f',pente_maxz),'�']);%    with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' �m] and \Lambda = ',sprintf('%3.2f',Lb),' �m']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
%saveas(gcf,['NULL_multi_MEAN_L.fig']);
print('-dpng', ['NULL_multi_MEAN_L.png'], '-r300');



temp=0;
for i=1:npts_pentez
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(i-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    load(['multi_slope_L_',sprintf('%3.2f',pente_deg),'.mat'])
    temp = temp+(nuldpt-mu_nu).^2;
end
s_nu = sqrt(temp./(npts_pentez-1));

% std
% ---
fmulti = s_nu;
dessin_multi
%set(gca,'xlim',[.3 .57])
set(gca,'xtick',[.3 .35 .4 .45 .5])
%set(gca,'ytick',[4 4.25 4.5 4.75 5 5.25])
%set(gca,'ylim',[4 5.25])
caxis([-3.8,-1]);  %large
set(hC,'ytick',[-3.5 -3 -2.5 -2 -1.5 -1])
%caxis([-4,-1]);  %zoom
%[vi,i]=min(fmulti,[],1);
%[min_fmulti,j]=min(vi);
%b=plot(Xparam(i(j)),Yparam(j),'Marker','+');
%set(b,'Linewidth',2,'color',[1 0 0]) 
t=title(['Standard Deviation (log_{10}) for \alpha = ',sprintf('%3.2f',pente_minz),'-',sprintf('%3.2f',pente_maxz),'�']);%     with B = [',sprintf('%3.1f',lb_min),' - ',sprintf('%3.1f',lb_max),' �m] and \Lambda = ',sprintf('%3.2f',Lb),' �m']);
set(t,'Fontname',fnz,'FontSize',fsz*1.2,'FontWeight',fwz,'HorizontalAlignment','center')
%saveas(gcf,['NULL_multi_STD_L.fig']);
print('-dpng', ['NULL_multi_STD_L.png'], '-r300');

