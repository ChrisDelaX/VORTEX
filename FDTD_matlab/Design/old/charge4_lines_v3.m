clear all
close all


lp = 4;
period = 1.4;%'1.42'
depth = 5.2;%'4.7'
fill = 0.45;%'0.41'
%swAngle = 3;%'6'
width = period*fill;%
%nWidths = 10;%'5'
nThetas = 16;%24;%64;%   ''multiple of 4 
nPeriods = 10;%30;%7;%
r2Max = period*nPeriods;
nLevels = 4;%8;%12;%

% Figures
% -------
newFig
grid off
set(gca,'visible','off')
axis equal

tic
depart_tmp=now;
depart=datestr(depart_tmp)


count=0;

nl = 1;
r2Min = 0;
rtmp = r2Max;
nThetastmp = nThetas;
nThetas = nThetas-2;
theta = 2*pi/nThetas;
widthtmp = width;
while nl <= nLevels    
%    for kk=0:ceil((nThetas-1)/4)
    for kk=0:(nThetas-1)
        nR = ceil(nPeriods*abs(cos(((lp/2-1)*kk-1/2*mySign(tan((lp/2-1)*kk*theta)))*theta))) ;        
        for ll=0:nR*nl-1
            width = widthtmp;
            draw = 1;
            r = width/2+ (ll)*period/abs(cos((lp/2-1)*kk*theta)) + mod((nl-1)*r2Max,period/abs(cos(((lp/2-1)*kk-1/2*mySign(tan((lp/2-1)*kk*theta)))*theta)));
            cosStart = cos(((lp/2-1)*kk+1/2)*theta);
            cosEnd = cos(((lp/2-1)*kk-1/2)*theta);
            rStart = r*cos((lp/2-1)*kk*theta)/cosStart;
            rEnd = r*cos((lp/2-1)*kk*theta)/cosEnd;
            if min(abs(rStart),abs(rEnd)) <= r2Max && max(abs(rStart),abs(rEnd)) >= r2Min
                xStart = rStart*cos((kk-1/2)*theta);
                yStart = rStart*sin((kk-1/2)*theta);
                xEnd = rEnd*cos((kk+1/2)*theta);
                yEnd = rEnd*sin((kk+1/2)*theta);
                startWedge = -((lp/2-1)*kk-1/2)*theta;
                endWedge = -((lp/2-1)*kk+1/2)*theta;
                if abs(rStart) >=  r2Max
                    if abs(rEnd) >=  r2Max - abs(width/2/cos(startWedge))
                        width = (r2Max-rEnd+abs(width/2/cos(startWedge)))*cos(startWedge);
                        rEnd = r2Max - abs(width/2/cos(startWedge));
                        draw = 2;
                    end
                    borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rEnd/r2Max*cosEnd))*mySign(sin((lp/2-1)*kk*theta));
                    xStart = r2Max*cos(borderAngle);
                    yStart = r2Max*sin(borderAngle);
                    endWedge = borderAngle-lp/2*kk*theta+pi/2;
                end
                if abs(rEnd) >=  r2Max
                    if abs(rStart) >=  r2Max - abs(width/2/cos(endWedge))
                        width = (r2Max-rStart+abs(width/2/cos(endWedge)))*cos(endWedge);
                        rStart = r2Max - abs(width/2/cos(endWedge));
                        draw = 2;
                    end
                    borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rStart/r2Max*cosStart))*mySign(sin((lp/2-1)*kk*theta));
                    xEnd = r2Max*cos(borderAngle);
                    yEnd = r2Max*sin(borderAngle);
                    startWedge = borderAngle-lp/2*kk*theta+pi/2;
                end
                if abs(rStart) <=  r2Min
                    if abs(rEnd) <=  r2Min + abs(width/2/cos(startWedge))
                        width = (rEnd-r2Min+abs(width/2/cos(startWedge)))*cos(startWedge);
                        rEnd = r2Min + abs(width/2/cos(startWedge));
                        draw = 3;
                    end
                    borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rEnd/r2Min*cosEnd))*mySign(sin((lp/2-1)*kk*theta));
                    xStart = r2Min*cos(borderAngle);
                    yStart = r2Min*sin(borderAngle);
                    endWedge = borderAngle-lp/2*kk*theta+pi/2;
                end
                if abs(rEnd) <=  r2Min
                    if abs(rStart) <=  r2Min + abs(width/2/cos(endWedge))
                        width = (rStart-r2Min+abs(width/2/cos(endWedge)))*cos(endWedge);
                        rStart = r2Min + abs(width/2/cos(endWedge));
                        draw = 3;
                    end
                    borderAngle = pi+(lp/2)*kk*theta+(pi/2+asin(rStart/r2Min*cosStart))*mySign(sin((lp/2-1)*kk*theta));
                    xEnd = r2Min*cos(borderAngle);
                    yEnd = r2Min*sin(borderAngle);
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
                if draw == 1  %grey
                    CreateTrapezoidWedged
                    count = count+1;
                elseif draw == 1%2  %red
                    CreateTrapezoidWedged2
                    count = count+1;
                elseif draw == 1%3  %blue
                    CreateTrapezoidWedged3
                    count = count+1;
                end
            end
            [xCell(1),yCell(1)] = pol2cart((kk+1/2)*theta,r2Max);
            [xCell(2),yCell(2)] = pol2cart((kk+1/2)*theta,r2Min);
            [xCell(3),yCell(3)] = pol2cart((kk-1/2)*theta,r2Min);
            [xCell(4),yCell(4)] = pol2cart((kk-1/2)*theta,r2Max);
            line(xCell,yCell,'color',mygreen)
            %patch(xCell,yCell,'k','facealpha',0,'edgecolor','k')
        end
    end
    
    nl = nl+1;
    r2Min = r2Min+rtmp;
    r2Max = r2Max+rtmp;
    nThetas = nThetas+nThetastmp;
    theta = 2*pi/nThetas;
    
end

count

r2Max = r2Max-rtmp;
%depart
%fin
toc
%datestr(now)

% Graphics
la_tmp=.05;%
%line([-r2Max r2Max]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
%line([-r2Max/sqrt(2) r2Max/sqrt(2)]*(1+la_tmp),[-r2Max/sqrt(2) r2Max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
%line([0 0],[-r2Max r2Max]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
%line([r2Max/sqrt(2) -r2Max/sqrt(2)]*(1+la_tmp),[-r2Max/sqrt(2) r2Max/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(r2Max*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-r2Max*(1+5*la_tmp),0,'$180^{\circ}$')
text(r2Max/sqrt(2)*(1+la_tmp),r2Max/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-r2Max/sqrt(2)*(1+4*la_tmp),-r2Max/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
text(-r2Max*la_tmp,r2Max*(1+2*la_tmp),'$90^{\circ}$')
text(-r2Max*la_tmp,-r2Max*(1+2*la_tmp),'$270^{\circ}$')
text(-r2Max/sqrt(2)*(1+4*la_tmp),r2Max/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(r2Max/sqrt(2)*(1+la_tmp),-r2Max/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
axis equal
xlabel('$\mu m$')
ylabel('$\mu m$')
set(gca,'xLim',[-r2Max r2Max]*(1+8*la_tmp))
set(gca,'ylim',[-r2Max r2Max]*(1+6*la_tmp))
%tick2latex

%save
print('-depsc2',sprintf('lines_sketch_lvl=%d_th=%d_r=%d.eps',nLevels,nThetastmp,nPeriods),'-r300');

% %Full screen
% pos0=get(0,'Screensize');
% pos1=get(gcf,'Position');
% pos2=[pos1(1)-pos1(3)*(pos0(4)/pos1(4)-1)/2 1 pos1(3)*pos0(4)/pos1(4) pos0(4)];
% set(gcf,'Position',pos2,'PaperUnits','points','PaperPosition',pos2);
% %set(gcf,'Position',get(0,'Screensize'),'PaperUnits','points','PaperPosition',get(0,'Screensize')); % Maximize figure + screen sized images [0 0 1440 878]



