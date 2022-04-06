close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 24; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width


% TRANS without filter
%--------------------


figure('name','trans')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'box','on')%,'linewidth',2)
set(gca,'XLim',[0.07 100])%,'xtick',bandtick)
set(gca,'XScale','log','XMinorGrid','on','XAxisLocation','top')
set(gca,'ylim',[0 26],'ytick',[],'ydir','reverse')
set(gca,'Fontname',fnz,'FontSize',0.8*fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
xlabel('Wavelength (\mum)','FontSize',0.8*fsz)
%ylabel('Transmission','FontSize',0.8*fsz)

mygrey = [.8 .8 .8];
% bands
longband=[.2 25.6;.2 25.6];
text(.08,1,'Bands:','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([0.39 .7],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(.35,1,'Visible','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([1.15,1.4],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(1.18,1,'J','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([1.5 1.8],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(1.5,1,'H','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([2 2.4],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(2,1,'K','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([3.5 4.1],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(3.45,1,'L''','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([4.5 5.1],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(4.45,1,'M','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
h=area([8 13],longband,'FaceColor',mygrey,'EdgeColor',mygrey,'linewidth',lwz*.01);
set(h(1),'visible','off')
text(9,1,'N','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
text(70,1,'FIR','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)


%Materials
go=[4 4];
step=2;
off=.7;
line([.200 3.4],go,'color','k','linewidth',lwz)
text(.200,go(1)-off,'SiO2','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.150 5],go,'color','k','linewidth',lwz)
text(.150,go(1)-off,'Al2O3','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.150 8],go,'color','k','linewidth',lwz)
text(.150,go(1)-off,'CaF2','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.150 8],go,'color','k','linewidth',lwz)
text(.150,go(1)-off,'MgF2','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.150 8],go,'color','k','linewidth',lwz)
text(.150,go(1)-off,'BaF2','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([1 8],go,'color','k','linewidth',lwz)
text(1,go(1)-off,'GaAs','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([1.06 6.7],go,'color','k','linewidth',lwz)
text(1.06,go(1)-off,'Si','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
line([30 100],go,'color','k','linewidth',lwz)
text(30,go(1)-off,'Si','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([2 17],go,'color','k','linewidth',lwz)
text(2,go(1)-off,'Ge','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.45 14],go,'color','k','linewidth',lwz)
text(.45,go(1)-off,'ZnS','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.5 20],go,'color','k','linewidth',lwz)
text(.5,go(1)-off,'ZnSe','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;
line([.1 100],go,'color','k','linewidth',lwz)
text(.11,go(1)-off,'CVD diamond','Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
go=go+step;




print('-depsc2',sprintf('TH_all_bands.eps'), '-r300');
