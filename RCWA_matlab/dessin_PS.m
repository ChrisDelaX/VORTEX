%lb_vect=[1.15 1.5 2 3.45 4.35 8.5 11];
lb_vect=[3.45 4.35 8.5 11];
npts_lb=size(lb_vect,2);
%npts_E=61;%11;%3;%
pente_deg = pente/pi*180;


close all
warning off MATLAB:polyfit

fnz = 'Arial'; % fontname
fsz = 20; % fontsize
fwz = 'normal';%'Bold'; % fontweight
msz = 8; % marker size
lwz = 2;  % line width
index = [1.5 2.8];
indextick = [1.5 1.75 2.0 2.25 2.5 2.75];
liss_index = [Xparam(1,1):1e-4:Xparam(1,end)];



figure('name','1: rms error')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % paper
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',index,'xtick',indextick)
%set(gca,'ylim',[0.4 1],'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Refractive index','FontSize',fsz*1.2)
ylabel('Phase-shift error (rms)','FontSize',fsz*1.2)
set(gca,'YScale','log')

for band_nb=1:npts_lb
    %lb=lb_vect(band_nb);
    PSerr_temp = smooth(log10(rms_err(band_nb,:)),100)';
    PSerr_temp(1) = PSerr_temp(2);
    PSerr_temp(end) = PSerr_temp(end-1);
    PSerr = 10.^polyval(polyfit(Xparam(1,:),PSerr_temp,9),liss_index);
    [b]=plot(liss_index,PSerr);%,'Marker','.');
    if band_nb == 1
        set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 2
        set(b,'color',[0 .4 .1],'linewidth',lwz)%,'marker','d','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 3
        set(b,'color',[.6 .1 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 4
        set(b,'color',[0 0 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    end
end
p=findobj('type','line');
idx=[4 3 2 1];
idy=[1 2 3 4];
str={' [ 3.5 - 4.1 \mum] L-band' ' [ 4.5 - 5.1 \mum] M-band' ' [ 8.0 - 9.6 \mum] lower N-band' ' [11.0 -13.2 \mum] upper N-band'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['rms_err.png'], '-r300');


%--------------------


figure('name','2: fill factor')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % Maximize figure + screen sized images
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',index,'xtick',indextick)
%set(gca,'ylim',[0.4 1],'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Refractive index','FontSize',fsz*1.2)
ylabel('Filling factor','FontSize',fsz*1.2)
%set(gca,'YScale','log')

for band_nb=1:npts_lb
    %lb=lb_vect(band_nb);
    [b]=plot(Xparam(1,:),smooth(Yparam1(band_nb,:),2));%,'Marker','.');
    if band_nb == 1
        set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 2
        set(b,'color',[0 .4 .1],'linewidth',lwz)%,'marker','d','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 3
        set(b,'color',[.6 .1 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 4
        set(b,'color',[0 0 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    end
end
p=findobj('type','line');
idx=[4 3 2 1];
idy=[1 2 3 4];
str={' [ 3.5 - 4.1 \mum] L-band' ' [ 4.5 - 5.1 \mum] M-band' ' [ 8.0 - 9.6 \mum] lower N-band' ' [11.0 -13.2 \mum] upper N-band'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['fill_fact.png'], '-r300');



%--------------------


figure('name','3: thickness')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % Maximize figure + screen sized images
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',index,'xtick',indextick)
%set(gca,'ylim',[0.4 1],'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Refractive index','FontSize',fsz*1.2)
ylabel('Grating depth (\mum)','FontSize',fsz*1.2)
%set(gca,'YScale','log')

for band_nb=1:npts_lb
    %lb=lb_vect(band_nb);
    thick_temp = smooth(Yparam3(band_nb,:),100)';
    thick_temp(1) = thick_temp(2);
    thick_temp(end) = thick_temp(end-1);
    thick = polyval(polyfit(Xparam(1,:),thick_temp,9),liss_index);
    [b]=plot(liss_index,thick);%,'Marker','.');
    if band_nb == 1
        set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 2
        set(b,'color',[0 .4 .1],'linewidth',lwz)%,'marker','d','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 3
        set(b,'color',[.6 .1 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 4
        set(b,'color',[0 0 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    end
end
p=findobj('type','line');
idx=[4 3 2 1];
idy=[1 2 3 4];
str={' [ 3.5 - 4.1 \mum] L-band' ' [ 4.5 - 5.1 \mum] M-band' ' [ 8.0 - 9.6 \mum] lower N-band' ' [11.0 -13.2 \mum] upper N-band'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['thickness.png'], '-r300');



%--------------------


figure('name','4: aspect ratio')
set(gcf,'color',[1 1 1])
%set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images
set(gcf,'Position',[400    0   800   800],'PaperUnits','points','PaperPosition',[400    0   800   800]); % Maximize figure + screen sized images
hold on
%grid on
set(gca,'Fontname',fnz,'FontSize',fsz,'FontWeight',fwz)
set(gca,'XLim',index,'xtick',indextick)
%set(gca,'ylim',[0.4 1],'ytick',[.4 .5 .6 .7 .8 .9 1])
set(gca,'ylim',[5 40])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')%,'YMinorGrid','on')
set(gca,'box','on','linewidth',2)
xlabel('Refractive index','FontSize',fsz*1.2)
ylabel('Aspect ratio','FontSize',fsz*1.2)
%set(gca,'YScale','log')

for band_nb=1:npts_lb
    %lb=lb_vect(band_nb);
    F_asp_rat = .5 - abs(.5-Yparam1(band_nb,:));
	asprat = Yparam3(band_nb,:)./(F_asp_rat.*Yparam2(band_nb,:));
    arliss_temp = smooth(asprat,8)';
    arliss_temp(1) = arliss_temp(2);
    arliss_temp(end) = arliss_temp(end-1);
    arliss = polyval(polyfit(Xparam(1,:),arliss_temp,12),liss_index);
    [b]=plot(liss_index,arliss);%,'Marker','.');
    %[b]=plot(Xparam(1,:),arliss_temp);%,'Marker','.');
    if band_nb == 1
        set(b,'color',[.8 .6 0],'linewidth',lwz)%,'marker','^','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 2
        set(b,'color',[0 .4 .1],'linewidth',lwz)%,'marker','d','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 3
        set(b,'color',[.6 .1 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    elseif band_nb == 4
        set(b,'color',[0 0 0],'linewidth',lwz)%,'marker','v','markersize',msz,'markerfacecolor','auto')
    end
end
p=findobj('type','line');
idx=[4 3 2 1];
idy=[1 2 3 4];
str={' [ 3.5 - 4.1 \mum] L-band' ' [ 4.5 - 5.1 \mum] M-band' ' [ 8.0 - 9.6 \mum] lower N-band' ' [11.0 -13.2 \mum] upper N-band'};
l=legend(p(idx),str(idy));
set(l,'Fontname',fnz,'FontSize',fsz*0.9,'location','northeast')
set(l,'box','on','linewidth',2)
print('-dpng', ['aspect_ratio.png'], '-r300');





