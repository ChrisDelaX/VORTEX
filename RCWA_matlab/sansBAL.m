f=figure('FileName','Diamant');
set(f,'color',[1 1 1])
hold on
%grid on

Fmoy = F+d/2*tan(pente)*2/Lb;
if Fmoy < 0.5
    F_aspRat = Fmoy;
else
    F_aspRat = 1-Fmoy;
end
aspRat = d/(F_aspRat*Lb);

nulldepth = subplot(2,2,[1 3]);
%plot(lb_t,nul_res_sp_b,'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
for i=1:size(nul_res_sp_b,2)
    null_res_sp(1,i)=mean(nul_res_sp_b);
end
meanRho=1/mean(nul_res_sp_b);
%plot(lb_t,polyval(polyfit(lb_t,nul_res_sp_b,9),lb_t),'k.-',lb_t,null_res_sp,'k--','Linewidth',2)
p1=plot(lb_t,nul_res_sp_b,'k.-',lb_t,null_res_sp,'k--','Linewidth',2);
set(p1,'Linewidth',2) 
axis square
set(nulldepth,'YScale','log')
xlabel('Wavelength \lambda (\mum)','FontSize',16,'FontWeight','bold')
title('AGPM Diamond : N-band','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')
l1=legend('Null Depth','mean Null Depth');
set(l1,'FontSize',16,'FontWeight','bold')
%gtext({['period \Lambda:  ',num2str(Lb),' \mum'],['total depth h:  ',num2str(d),' \mum'],['AR thickness h_{AR}:  ',num2str(d_AR),' \mum'],['filling factor F:  ',num2str(F*100),' %'],['aspect ratio A-R:  ',num2str(aspRat)],['mean null depth µ:  ',num2str(mean(nul_res_sp_b))]})
gtext({['period \Lambda:  ',num2str(Lb),' \mum'],['depth h:  ',num2str(d),' \mum'],['slope \alpha:  ',num2str(pente/pi*180),' °'],['filling factor F:  ',num2str(F*100),' %'],['aspect ratio A-R:  ',num2str(aspRat)],['mean null depth µ:  ',num2str(mean(nul_res_sp_b))]},'FontSize',16,'FontWeight','bold')
%gtext({['period \Lambda:  ',num2str(Lb),' \mum'],['depth h:  ',num2str(d),' \mum'],['mean null depth µ:  ',num2str(mean(nul_res_sp_b))]})



% dessinprofil = subplot(2,3,2);
% imagesc(prof)
% xlabel('Period \Lambda (%)')
% title('Grating thickness (%)')

% limiteZOG = subplot(2,3,3);
% plot(lb_t,limZOG,lb_t,Lb,'k-.','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% ylabel('Limite ZOG \lambda/n(\lambda) (\mum)')

transmittance = subplot(2,2,4);
tmp12=mean(tmp1+tmp2)/2;
%TsansZOG=2./(n2lb+1);
%plot(lb_t,tmp1,'-',lb_t,tmp2,'--',lb_t,tmp12,'k:',lb_t,TsansZOG,'k-.','Linewidth',2)
p2=plot(lb_t,tmp1,'-',lb_t,tmp2,'--',lb_t,tmp12,'k:','Linewidth',2);
set(p2,'Linewidth',2) 
xlabel('Wavelength \lambda (\mum)','FontSize',16,'FontWeight','bold')
title('Transmittance','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')
l2=legend('TE','TM');
set(l2,'FontSize',16,'FontWeight','bold')

% reflectance = subplot(2,2,4);
% tmp45=mean(tmp4+tmp5)/2;
% %RsansZOG=(n2lb-1)./(n2lb+1);
% %plot(lb_t,tmp4,'-',lb_t,tmp5,'--',lb_t,tmp45,'k:',lb_t,RsansZOG,'k-.','Linewidth',2)
% plot(lb_t,tmp4,'-',lb_t,tmp5,'--',lb_t,tmp45,'k:','Linewidth',2)
% xlabel('Wavelength \lambda (\mum)')
% title('Reflectance')
% legend('TE','TM')

phaseshift = subplot(2,2,2);
for i=1:size(tmp3,2)
    tmp3moy(1,i)=mean(tmp3);
    tmpPi(1,i)=pi;
end
p3=plot(lb_t,tmp3,lb_t,tmp3moy,'k:',lb_t,tmpPi,'k-.','Linewidth',2);
set(p3,'Linewidth',2) 
xlabel('Wavelength \lambda (\mum)','FontSize',16,'FontWeight','bold')
title('TE-TM Phase Shift  \Delta\Phi_{TE-TM}  (rad)','FontSize',16,'FontWeight','bold')
set(gca,'FontSize',16,'FontWeight','bold')
l3=legend('\Delta\Phi_{TE-TM}','mean \Delta\Phi_{TE-TM}','\pi = 3.14 rad');
set(l3,'FontSize',16,'FontWeight','bold')

