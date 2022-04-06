pente_minz = 2.6;
pente_maxz = 2.8;
npts_pentez = 11;%5;%3;%

F_minz=0.35;
F_maxz=0.45;
npts_Fz=3;
% 
% for k=1:npts_pentez
%     
%     pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
%     pente=rad(pente_deg);
%     
%     for i=1:npts_Fz
%         
%         F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);
%         
%         load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
%         dessin_spectrum_RMS
%         
%     end
% end
% 

close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 20; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width


band = [11 13.2];
bandtick = [11.0 11.5 12.0 12.5 13.0];
liss_lb_t = [lb_t(1):1e-4:lb_t(end)];




figure('name','transmittance')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'ylim',[80 100],'ytick',[80 85 90 95 100])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('Transmission (%)','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);

        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
%            b=plot(lb_t,100.*(tmp1+tmp2)./2,'color','k');
            if pente==rad(2.6)
                %set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
            b=plot(lb_t,100.*(tmp1+tmp2)./2,'color','k','linewidth',lwz*1.5);
                %set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.80)
                %set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
TsansZOG= n2lb.*(2./(n2lb+1)).^2;
plot(lb_t,100.*TsansZOG,'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[12 11 7 6 2 1];
idy=[1 2 3 4 5 6];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°' ' without zog'};
%l=legend(p(idx),str(idy));
l=legend(' alpha ~2.6-2.8°',' without zog');
set(l,'Fontname',fnz,'FontSize',fsz*0.9);%,'location','southeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['Trans.png'], '-r300');




figure('name','TETMratio')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'ylim',[0.9 1.02]);%,'ytick',[.9 .95 1 1.05 1.1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('TE/TM ratio','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);

        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
            b=plot(lb_t,(tmp1./tmp2),'color','k');
            if pente==rad(2.6)
                set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
                set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.8)
                set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
plot(lb_t,zeros(size(lb_t)) + 1,'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[12 11 7 6 2 1];
idy=[1 2 3 4 5 6];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°' ' TE = TM'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['TETMratio.png'], '-r300');




figure('name','TETMdiff')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',band,'xtick',bandtick)
%set(gca,'YScale','log')
set(gca,'ylim',[-.1 .02])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('TE-TM differential','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);

        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
            b=plot(lb_t,tmp1-tmp2,'color','k');
            if pente==rad(2.6)
                set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
                set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.8)
                set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
plot(lb_t,zeros(size(lb_t)),'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[12 11 7 6 2 1];
idy=[1 2 3 4 5 6];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°' ' TE = TM'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['TETMdiff.png'], '-r300');



figure('name','logTETMdiff')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'YScale','log')
set(gca,'ylim',[.01 1])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('TE-TM differential (log_{10})','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);

        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
            TETMdiff = polyval(polyfit(lb_t,tmp1-tmp2,8),liss_lb_t);
            b=plot(liss_lb_t,abs(TETMdiff),'color','k');
            if pente==rad(2.6)
                set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
                set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.8)
                set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
%plot(lb_t,zeros(size(lb_t)) + 1,'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[11 10 6 5 1];
idy=[1 2 3 4 5];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','best')
set(l,'box','on','linewidth',2)
print('-dpng', ['logTETMdiff.png'], '-r300');





figure('name','retardance')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);
        
        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
            [AX,H1,H2] = plotyy(lb_t,retard,lb_t,retard*360);
            %retmoy = zeros(size(lb_t)) + mean(retard);
            %plot(lb_t,retmoy,'k--','Linewidth',1)
            %axis tight
            set(AX(1),'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz,'XLim',band,'xtick',bandtick,'XMinorTick','on','YColor','k','ylim',[0.47 0.53],'ytick',[])
            set(AX(2),'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz,'XLim',band,'xtick',bandtick,'XMinorTick','on','YColor','k','ylim',get(AX(1),'ylim').*360,'ytick',[])
            set(H1,'color','k')
            set(H2,'color','k')
            if pente==rad(2.6)
                set(H1,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
                set(H2,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
                set(H1,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
                set(H2,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.8)
                set(H1,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
                set(H2,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
plot(lb_t,zeros(size(lb_t)) + .5,'k-.','linewidth',lwz)
ylabel(AX(1),'Retardance (Waves)','FontSize',fsz*1.2)
ylabel(AX(2),'Phase Shift (°)','FontSize',fsz*1.2)
set(AX(1),'ylim',[0.47 0.53],'ytick',[.47 .48 .49 .50 .51 .52 .53])
set(AX(2),'ylim',get(AX(1),'ylim').*360,'ytick',[170 175 180 185 190])
p=findobj('type','line');
idx=[11 10 6 5 1 12];
idy=[1 2 3 4 5 6];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°' ' ideal HWP'};
l=legend(p(idx),str(idy));
%l=legend(' \alpha = 2.60°', '  \alpha = 2.70°', ' \alpha = 2.80°', 'ideal HWP','location','northeast');
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['retardance.png'], '-r300');





figure('name','logretard')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878] 
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper 
hold on
%grid on
set(gca,'XLim',band,'xtick',bandtick)
set(gca,'YScale','log')
set(gca,'ylim',[.0001 .2])%,'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Wavelength (\mum)','FontSize',fsz*1.2)
ylabel('| \Delta\Phi_{TE-TM} - \pi|  (log_{10})','FontSize',fsz*1.2)

for k=1:npts_pentez
    
    pente_deg = pente_minz+(pente_maxz-pente_minz)*(k-1)/(npts_pentez-1);
    pente=rad(pente_deg);
    
    for i=1:npts_Fz
        F=F_minz+(F_maxz-F_minz)*(i-1)/(npts_Fz-1);

        if F==0.4
            load(['spectrum_',sprintf('%3.2f',pente_deg),'_F_',sprintf('%3.2f',F),'.mat'])
            logretard = polyval(polyfit(lb_t,(retard*2*pi-pi),8),liss_lb_t);
            b=plot(liss_lb_t,abs(logretard),'color','k');
            if pente==rad(2.6)
                set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.70)
                set(b,'color',[.0 .4 .1],'linewidth',lwz*2)%,'marker','d','markersize',msz,'markerfacecolor','auto')
            elseif pente==rad(2.8)
                set(b,'color',[.6 .1 .0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
            end
        end
    end
end
%plot(lb_t,zeros(size(lb_t)) + 1,'k-.','linewidth',lwz)
p=findobj('type','line');
idx=[11 10 6 5 1];
idy=[1 2 3 4 5];
str={' alpha = 2.60°' ' ...' ' alpha = 2.70°' ' ...' ' alpha = 2.80°'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','east')
set(l,'box','on','linewidth',2)
print('-dpng', ['logretard.png'], '-r300');


