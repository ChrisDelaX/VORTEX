figure('FileName','Diamant')
hold on
%grid on
nulldepth = subplot(2,4,1:2);
plot(lb_t,nul_res_sp_b,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
axis square
set(nulldepth,'YScale','log')
xlabel('Wavelength \lambda (\mum)')
ylabel('Null Depth')
title('AGPM Diamond : N-band')
gtext({['Slope \alpha: ',num2str(pente/pi*180),' °'],['Period \Lambda: ',num2str(Lb),' \mum'],['Total Grating thickness d: ',num2str(d),' \mum'],['Filling factor F: ',num2str(F*100),' %']})
transmittance = subplot(2,4,3);
tmp12=mean(tmp1+tmp2)/2;
TsansZOG=2./(n2lb+1);
plot(lb_t,tmp1,'-',lb_t,tmp2,'--',lb_t,tmp12,'k:',lb_t,TsansZOG,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('Transmittance')
legend('TE','TM')
reflectance = subplot(2,4,4);
tmp45=mean(tmp4+tmp5)/2;
RsansZOG=(n2lb-1)./(n2lb+1);
plot(lb_t,tmp4,'-',lb_t,tmp5,'--',lb_t,tmp45,'k:',lb_t,RsansZOG,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('Reflectance')
legend('TE','TM')

dessinprofil = subplot(2,4,5);
imagesc(prof)
xlabel('Period \Lambda (%)')
ylabel('Grating thickness (%)')

phaseshift = subplot(2,4,6);
tmp3moy=mean(tmp3);
plot(lb_t,tmp3,lb_t,tmp3moy,'k:',lb_t,pi,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('TE-TM Phase Shift  \Delta\Phi_{TE-TM}  (rad)')
limiteZOG = subplot(2,4,7);
plot(lb_t,limZOG,lb_t,Lb,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('Limite ZOG \lambda/n(\lambda) (\mum)')
autresordres = subplot(2,4,8);
nonZOG = 1-(tmp1+tmp2)/2-(tmp4+tmp5)/2;
plot(lb_t,nonZOG,'k-.','Linewidth',2)
xlabel('Wavelength \lambda (\mum)')
ylabel('Ordres non nuls')

