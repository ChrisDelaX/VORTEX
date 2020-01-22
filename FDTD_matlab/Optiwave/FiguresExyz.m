% Ex Real
% -------
newFig
hS = surfc(realGridFig,realGridFig',real(rot90(ExFig,2)));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\Re(\bf{E_x})$')
tick2latex
%colorbar2latex(hC)
print('-depsc2',sprintf('%s/ExReal.eps',sgvc), '-r300')

% Ex Imag
% -------
newFig
hS = surfc(realGridFig,realGridFig',imag(rot90(ExFig,2)));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\Im(\bf{E_x})$')
tick2latex
%colorbar2latex(hC)
print('-depsc2',sprintf('%s/ExImag.eps',sgvc), '-r300')

% Ey Real
% -------
newFig
hS = surfc(realGridFig,realGridFig',real(rot90(EyFig,2)));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\Re(\bf{E_y})$')
tick2latex
%colorbar2latex(hC)
print('-depsc2',sprintf('%s/EyReal.eps',sgvc), '-r300')

% Ey Imag
% -------
newFig
hS = surfc(realGridFig,realGridFig',imag(rot90(EyFig,2)));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz)
set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('$\Im(\bf{E_y})$')
tick2latex
%colorbar2latex(hC)
print('-depsc2',sprintf('%s/EyImag.eps',sgvc), '-r300')
% 
% % Ez Real
% % -------
% newFig
% hS = surfc(realGridFig,realGridFig',real(rot90(EzFig,2)));
% set(hS,'EdgeColor','none')
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz)
% set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$\Re(\bf{E_z})$')
% tick2latex
% %colorbar2latex(hC)
% print('-depsc2',sprintf('%s/EzReal.eps',sgvc), '-r300')
% 
% % Ez Imag
% % -------
% newFig
% hS = surfc(realGridFig,realGridFig',imag(rot90(EzFig,2)));
% set(hS,'EdgeColor','none')
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz)
% set(hC,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$\Im(\bf{E_z})$')
% tick2latex
% %colorbar2latex(hC)
% print('-depsc2',sprintf('%s/EzImag.eps',sgvc), '-r300')