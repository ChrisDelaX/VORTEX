q=(tmp1./tmp2);

nul_phase_fin=(tmp3-pi).^2;

%resultat du nulling avec phase et amplitude mismatch
%----------------------------------------------------
% nul_res_sp_b(ib)=nul_phase_fin(ib)/4;
nul_res_sp_b=(((1+sqrt(q)).^2)./((1-sqrt(q)).^2+nul_phase_fin.*sqrt(q))).^(-1);



figure;
subplot(2,2,1);
plot(lb_t,nul_res_sp_b,'color',[0 0 0],'Linewidth',2);hold on;plot(lb_t,tmp1(:).*0.005.*tmp4(:).*0.995+tmp2(:).*0.005.*tmp5(:).*0.995,'color',[0 0 0],'Linewidth',2,'Linestyle',':');
xlabel('Wavelength (microns)');ylabel('Null depth/ghost');axis square;set(gca,'YScale','log');

subplot(2,2,2);plot(lb_t,tmp1,lb_t,tmp2);axis square;xlabel('Wavelength (microns)');ylabel('Transmittance');
subplot(2,2,3);plot(lb_t,tmp3);axis square;xlabel('Wavelength (microns)');ylabel('Phase shift');
subplot(2,2,4);imagesc(prof);

% figure(1);set(1,'position',[256 256 512 488]);
% % subplot(1,2,1);
% plot(lb_t,tmp3);xlabel('Wavelength (microns)');ylabel('Phase shift');axis square;
% % subplot(1,2,2);
% figure(2);set(2,'position',[256 256 512 488]);
% plot(lb_t,tmp1,lb_t,tmp2);xlabel('Wavelength (microns)');ylabel('Transmittance');axis square;


