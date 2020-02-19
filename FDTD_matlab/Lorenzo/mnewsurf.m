%% based on newsurf.m by Crhsitan
% sets up general setting for plotting the structure or the E field
% needs xMesh, yMesh and rez, vis (show the plot or not)

h=figure('visible',vis);
hold on
grid on
mywhite = [1 1 1];
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
set(colorbar,'box','on','linewidth',2)
axis equal
%caxis([0 12]) %specify caxis range, generally not necessary
axis([1 xMesh 1 yMesh])
xlabel(['x-mesh = ' num2str(1/rez) ' \mum (' num2str(xMesh) 'pts)'])
ylabel(['y-mesh = ' num2str(1/rez) ' \mum (' num2str(yMesh) 'pts)'])
