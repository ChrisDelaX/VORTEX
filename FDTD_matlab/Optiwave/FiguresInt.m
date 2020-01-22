% % Ax
% % --
% newFig
% hS = surfc(realGridFig,realGridFig',rot90(AxFig,2));
% set(hS,'EdgeColor','none')
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz,'FontSize',fsz) 
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$A_x=|\bf{E_x}\:|$')
% tick2latex
% %colorbar2latex(hC)
% print('-depsc2',sprintf('%s/Ax.eps',sgvc), '-r300')
% 
% % Ay
% % --
% newFig
% hS = surfc(realGridFig,realGridFig',rot90(AyFig,2));
% set(hS,'EdgeColor','none')
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz,'FontSize',fsz) 
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$A_y=|\bf{E_y}\:|$')
% tick2latex
% %colorbar2latex(hC)
% print('-depsc2',sprintf('%s/Ay.eps',sgvc), '-r300')
% 
% % Az
% % --
% newFig
% hS = surfc(realGridFig,realGridFig',rot90(AzFig,2));
% set(hS,'EdgeColor','none')
% %shading('interp')   %  ONLY if meshNumber < 70 !!!
% hC = colorbar;
% set(hC,'box','on','linewidth',lwz,'FontSize',fsz) 
% %graphics
% axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
% xlabel('X axis $(\mu m)$')
% ylabel('Y axis $(\mu m)$')
% title('$A_z=|\bf{E_z}\:|$')
% tick2latex
% %colorbar2latex(hC)
% print('-depsc2',sprintf('%s/Az.eps',sgvc), '-r300')

% Intensity I
% -----------
newFig
%fftInt=abs(otf2psf(Int));
%plot(realGrid,fftInt(:,226))
hS = surfc(realGridFig,realGridFig',rot90(IntFig,2));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('Intensity $I=A^2=A^2_x+A^2_y+A^2_z$')
tick2latex
print('-depsc2',sprintf('%s/Int.eps',sgvc), '-r300')


% Intensity I log scale
% ---------------------
newFig
hS = surfc(realGridFig,realGridFig',rot90(log10(IntFig),2));
set(hS,'EdgeColor','none')
%shading('interp')   %  ONLY if meshNumber < 70 !!!
hC = colorbar;
set(hC,'box','on','linewidth',lwz,'FontSize',fsz)
%set(get(hC,'YLabel'),'String','log$_{10}$','Rotation',-90,'FontSize',fsz)
%graphics
axis([realGridFig(1) realGridFig(end) realGridFig(1) realGridFig(end)])
xlabel('X axis $(\mu m)$')
ylabel('Y axis $(\mu m)$')
title('Intensity log$_{10}(I)$')
tick2latex
colorbar2latex(hC)
print('-depsc2',sprintf('%s/Intlog.eps',sgvc), '-r300')


% mean intensity
% --------------
newFig
grid on
plot(offAxis,meanInt,'-','color',mygreen,'linewidth',lwz)
%graphics
xlabel('Off-axis distance $(\mu m)$')
ylabel('Intensity $I$')
title('Mean radial (off-axis) intensity')
%leg=legend(' temp')
%set(leg,'box','on','location','northwest')
tick2latex
print('-depsc2',sprintf('%s/meanInt.eps',sgvc), '-r300')

