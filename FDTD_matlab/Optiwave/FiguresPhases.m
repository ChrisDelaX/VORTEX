% phix
% ----
newFig
hS = surfc(realGridFig,realGridFig',rot90(phixFig./pi,2));
set(hS,'EdgeColor','none')
caxis([0 2])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\phi_x = \mathrm{arg}\:\:\bf{E_x} \: \: (\pi rad)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/phix.eps',sgvc), '-r300')

% phiy
% ----
newFig
hS = surfc(realGridFig,realGridFig',rot90(phiyFig./pi,2));
set(hS,'EdgeColor','none')
caxis([0 2])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz) 
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\phi_y = \mathrm{arg}\:\:\bf{E_y} \: \: (\pi rad)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/phiy.eps',sgvc), '-r300')
% 
% % phiz
% % ----
% newFig
% hS = surfc(realGridFig,realGridFig',rot90(phizFig./pi,2));
% set(hS,'EdgeColor','none')
% caxis([0 2])
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz,'FontSize',fsz) 
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$\phi_z = \mathrm{arg}\:\:\bf{E_z} \: \: (\pi rad)$')
% tick2latex
% colorbar2latex(hC)
% print('-depsc2',sprintf('%s/phiz.eps',sgvc), '-r300')

% phitetm
% -------
newFig
hS = surfc(realGridFig,realGridFig',rot90(phitetmFig./pi,2));
set(hS,'EdgeColor','none')
caxis([0 2])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\phi_{TE-TM} \:\: = \mathrm{arg}\:\:\bf{E_{TM}} \:\: - \mathrm{arg}\:\:\bf{E_{TE}} \: \: (\pi rad)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/phiTETM.eps',sgvc), '-r300')

% Residuals phitetm
% -----------------
newFig
hS = surfc(realGridFig,realGridFig',rot90(abs(phitetmFig-pi)./pi,2));
set(hS,'EdgeColor','none')
caxis([0 1])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\phi_{TE-TM} \:$   Residuals $(\pi rad)$')
tick2latex
%colorbar2latex(hC)
print('-depsc2',sprintf('%s/residPhiTETM.eps',sgvc), '-r300')

% Pancharatnam
% ------------
newFig
hS = surfc(realGridFig,realGridFig',rot90(phiPanFig./pi,2));
set(hS,'EdgeColor','none')
caxis([0 2])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('Pancharatnam phase $\phi_p \: \: (\pi rad)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/phiPan.eps',sgvc), '-r300')

% Residuals Pancha  
% ----------------
newFig
residPAN = mod(phiPanFig+orientFig.*2+pi,2*pi)-pi;
hS = surfc(realGridFig,realGridFig',rot90(abs(residPAN)./pi,2));
set(hS,'EdgeColor','none')
caxis([0 1])
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('Pancharatnam Residuals $(\pi rad)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/residPhiPan.eps',sgvc), '-r300')


