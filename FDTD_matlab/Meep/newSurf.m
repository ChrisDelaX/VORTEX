% by Christian Delacroix (2014)

%close all
figure
hold on
grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.4,0.4,0.4];
fnz = 'Arial'; % fontname: Helvetica
fsz = 12; % fontsize: 10
fwz = 'normal';  % fontweight: bold
msz = 8; % marker size
lwz = 2;  % line width
set(gcf,'color',mywhite)
set(gca,'box','on','linewidth',lwz)
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
%set(gca,'ticklength',-.9*get(gca,'ticklength'))
hC = colorbar;
set(hC,'box','on','linewidth',2)
axis equal
xlabel(['x-mesh = ' num2str(1/rez) ' \mum (' num2str(xMesh) 'pts)'])
ylabel(['y-mesh = ' num2str(1/rez) ' \mum (' num2str(yMesh) 'pts)'])

%set(0,'DefaultTextInterpreter', 'latex')



% % Vertically stretch to screen size
% pos0=get(0,'Screensize');
% pos1=get(gcf,'Position');
% pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
% set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);

% % Maximize figure + screen sized images [0 0 1440 878]
% set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); 
