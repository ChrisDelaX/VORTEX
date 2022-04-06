clear all
close all


lp = 4;
period = 1.4;%'1.42'
depth = 5.2;%'4.7'
fill = 0.45;%'0.41'
%swAngle = 3;%'6'
width = period*fill;%
%nWidths = 10;%'5'
nThetas = 62;%14;%   ''multiple of 4 minus 2
theta = 2*pi/nThetas;
nPeriods = 28;%7;%
rMax = period*nPeriods;

 
% Figures
% -------
figure
hold on
%grid on
mywhite = [1 1 1];
mygreen = [0 .5 0];
myred = [1 .2 0];
myblue = [0 .2 1];
mygrey = [0.35,0.35,0.35];
%mygrey = [.55,.55,.55];
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
set(0,'DefaultTextInterpreter', 'latex')
set(gca,'visible','off')
axis equal

tic
depart_tmp=now;
depart=datestr(depart_tmp)
    



for kk=0:nThetas-1
    nR = ceil(nPeriods*abs(cos(((lp/2-1)*kk-1/2*mySign(tan((lp/2-1)*kk*theta)))*theta)));
    for ll=0:nR-1
        r = (ll+1/2)*period/abs(cos((lp/2-1)*kk*theta));
        cosStart = cos(((lp/2-1)*kk+1/2)*theta);
        cosEnd = cos(((lp/2-1)*kk-1/2)*theta);
        rStart = r*cos((lp/2-1)*kk*theta)/cosStart;
        rEnd = r*cos((lp/2-1)*kk*theta)/cosEnd;
        if min(abs(rStart),abs(rEnd)) < rMax
            xStart = rStart*cos((kk-1/2)*theta);
            yStart = rStart*sin((kk-1/2)*theta);
            xEnd = rEnd*cos((kk+1/2)*theta);
            yEnd = rEnd*sin((kk+1/2)*theta);
            startWedge = -((lp/2-1)*kk-1/2)*theta;
            endWedge = -((lp/2-1)*kk+1/2)*theta;
            if abs(rStart) >=  rMax
                borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rEnd/rMax*cosEnd))*mySign(sin((lp/2-1)*kk*theta));
                xStart = rMax*cos(borderAngle);
                yStart = rMax*sin(borderAngle);
                endWedge = borderAngle-lp/2*kk*theta+pi/2;
            end
            if abs(rEnd) >=  rMax
                borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rStart/rMax*cosStart))*mySign(sin((lp/2-1)*kk*theta));
                xEnd = rMax*cos(borderAngle);
                yEnd = rMax*sin(borderAngle);
                startWedge = borderAngle-lp/2*kk*theta+pi/2;
            end
            if mySign(cos(((lp/2)-1)*kk*theta)) == -1
                temp = startWedge;
                startWedge = endWedge;
                endWedge = temp;
            end
            length = sqrt((xStart-xEnd)^2+(yStart-yEnd)^2);
            xPos = (xStart+xEnd)/2;
            yPos = (yStart+yEnd)/2;
            twist = kk*theta*lp/2;
            CreateTrapezoidWedged
        end
    end
end



%depart
%fin
toc
%datestr(now)

% Graphics
la_tmp=.05;%
line([-rMax rMax]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
text(rMax*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-rMax*(1+5*la_tmp),0,'$180^{\circ}$')
line([-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(rMax/sqrt(2)*(1+la_tmp),rMax/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),-rMax/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
line([0 0],[-rMax rMax]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-rMax*la_tmp,rMax*(1+2*la_tmp),'$90^{\circ}$')
text(-rMax*la_tmp,-rMax*(1+2*la_tmp),'$270^{\circ}$')
line([rMax/sqrt(2) -rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-rMax/sqrt(2)*(1+4*la_tmp),rMax/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),-rMax/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
axis equal
xlabel('$\mu m$')
ylabel('$\mu m$')
set(gca,'xLim',[-rMax rMax]*(1+8*la_tmp))
set(gca,'ylim',[-rMax rMax]*(1+6*la_tmp))
%tick2latex

%save
print('-depsc2',sprintf('lines_sketch_th=%d_r=%d.eps',nThetas,nPeriods), '-r300');

% %Full screen
% pos0=get(0,'Screensize');
% pos1=get(gcf,'Position');
% pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
% set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



