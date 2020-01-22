clear all
close all

% INPUT PARAMETERS
% ----------------
periodStart = 1.42;%1;%
periodEnd = 1.13;%2*depth*tan(swa)/(1-fillEnd);
fillStart = 0.46; %0.2;%
fillEnd = 0.49; %0.5;%
lp = 4;     % topological charge
nbQuadrants = 1;  % 1,2 or 4
drawCurves = 1;  % 0 or 1
drawLines = 0;  % 0 or 1
drawCells = 0;  % 0 or 1
thDisc = deg2rad(5); % curves discretization (in rad) --> nThCurves
nThLines = 7;%10;%  % lines discretization (number)
offset = 5;%2;%  % between 2 levels (in periods)
rMax = 100;%5;%5000;%

% merging verification
% --------------------
depth = 5.3;% %grating depth
swa = deg2rad(0);   % sidewall angle
if periodEnd < 2*depth*tan(swa)/(1-fillEnd)
    2*depth*tan(swa)/(1-fillEnd)
    error('periodRat too small. Merging of the walls.')
end

% Initializing figure and values
% ------------------------------
tic
datestr(now)
newFig
grid off
set(gca,'visible','off')
axis equal
periodRat = periodEnd/periodStart;
fillRat = fillEnd/fillStart;
cosStart = 1;
cosEnd = periodRat^(2/lp);
thStart = acos(cosStart);
thEnd = acos(cosEnd);
thMoy = (thStart+thEnd)/2;
stepPhi = periodStart;
stepThStart = stepPhi*cos(thStart);
stepThEnd = stepPhi*cos(thEnd);
rMin = 0;
jj=0;
outOfBounds = 0;
while outOfBounds==0
    
    % Curves along x-axis
    % -------------------
    nbStart = ceil(rMin/cos((thEnd-thStart)/2)/stepThStart+1/2);
    nbEnd = ceil(rMax/stepThEnd-1/2);
    for ii=nbStart:nbEnd
        % Draws Cells
        if ii==nbStart
            clear xCell yCell
            rCell=rMin/cos((thEnd-thStart)/2);
            [xCell(1),yCell(1)] = pol2cart(thStart,rMax);
            [xCell(2),yCell(2)] = pol2cart(thStart,rCell);
            [xCell(3),yCell(3)] = pol2cart(thEnd,rCell);
            [xCell(4),yCell(4)] = pol2cart(thEnd,rMax);
            if drawCells == 1, line(xCell,yCell,'color',myblue,'linewidth',lwz)
                if nbQuadrants >= 2, line(-xCell,yCell,'color',myblue,'linewidth',lwz)
                end
                if nbQuadrants == 4,
                    line(-xCell,-yCell,'color',myblue,'linewidth',lwz)
                    line(xCell,-yCell,'color',myblue,'linewidth',lwz)
                end
            end
        end
        % Draws Curves
        rPhi = stepPhi*(2*ii-1)/4;
        rThStart = 2*rPhi*cos(thStart);
        rThEnd = 2*rPhi*cos(thEnd);
        phiStart0=2*thStart;
        phiEnd0=2*thEnd;
        nThCurves = ceil((phiEnd0-phiStart0)/thDisc);
        for kk=1:nThCurves
            phiStart = phiStart0+(phiEnd0-phiStart0)/nThCurves*(kk-1);
            phiEnd = phiStart0+(phiEnd0-phiStart0)/nThCurves*(kk);
            widthStart = fillStart*periodStart;
            widthEnd = fillEnd*periodEnd;
            width = widthStart+(widthEnd-widthStart)/nThCurves*(kk-1/2);
            twist = (phiStart+phiEnd)/2;
            length = rPhi*(phiEnd-phiStart);
            [xPos,yPos] = pol2cart(twist,rPhi);
            xPos = xPos+rPhi;
            [thPos,rPos] = cart2pol(xPos,yPos);
            startWedge = thStart+(thEnd-thStart)/nThCurves*(kk-1)-twist;
            endWedge = thStart+(thEnd-thStart)/nThCurves*(kk)-twist;
            if rPos<=rMax && rPos>=rMin && drawCurves == 1, PatchWedgedBlock
                if nbQuadrants >= 2, PatchWedgedBlockQ2
                end
                if nbQuadrants == 4, PatchWedgedBlockQ4
                end
            end
        end
    end
    
    
    % increments
    % ----------
    jj=jj+1;
    cosStart = periodRat^(2*jj/lp);
    cosEnd = periodRat^(2*(jj+1)/lp);
    thStart = acos(cosStart);
    thStart2 = pi/4-abs(thStart-pi/4);
    thEnd = acos(cosEnd);
    thEnd2 = pi/4-abs(thEnd-pi/4);
    thMoy = (thStart+thEnd)/2;
    stepPhi = periodStart/periodRat^(jj);
    stepThStart = stepPhi*cos(thStart);
    stepThEnd = stepPhi*cos(thEnd);
    r2Min = rMin;
    rMin = offset/2/tan((thEnd-thStart)/2);
    r2Max = min(rMin,rMax);
    
    % lines along y-axis
    % ------------------
    theta = (pi/2-thStart)/nThLines;
    nPeriods = ceil(r2Max/stepThStart+1/2);
    
    for kk=1:nThLines
        % Draws Cells
        th2Start = thStart+(pi/2-thStart)/nThLines*(kk-1);
        th2End = thStart+(pi/2-thStart)/nThLines*(kk);
        th2Moy = (th2Start+th2End)/2;
        [xCell(1),yCell(1)] = pol2cart(th2End,r2Max);
        [xCell(2),yCell(2)] = pol2cart(th2End,r2Min);
        [xCell(3),yCell(3)] = pol2cart(th2Start,r2Min);
        [xCell(4),yCell(4)] = pol2cart(th2Start,r2Max);
        if drawCells == 1, line(xCell,yCell,'color',mygreen,'linewidth',lwz)
            if nbQuadrants >= 2, line(-xCell,yCell,'color',mygreen,'linewidth',lwz)
            end
            if nbQuadrants == 4,
                line(-xCell,-yCell,'color',mygreen,'linewidth',lwz)
                line(xCell,-yCell,'color',mygreen,'linewidth',lwz)
            end
        end
        % Draws Lines
        twist = lp/2*th2Moy;
        ll=0;
        for ll=0:nPeriods
            width = widthStart;
            r = (ll+1/2)*periodStart/cos(th2Moy) + mod((jj-1)*r2Max,periodStart/cos(th2Start));
            cosStart = cos(th2End);
            cosEnd = cos(th2Start);
            rStart = r*cos(th2Moy)/cosStart;
            rEnd = r*cos(th2Moy)/cosEnd;
            if min(rStart,rEnd) <= r2Max && max(rStart,rEnd) >= r2Min
                [xStart,yStart] = pol2cart(th2Start,rStart);
                [xEnd,yEnd] = pol2cart(th2End,rEnd);
                startWedge = th2Start-twist;
                endWedge = th2End-twist;
                if rStart >= r2Max
                    if rEnd >= r2Max-width/2/cos(endWedge)
                        width = (r2Max-rEnd)*cos(endWedge)+width/2;
                        rEnd = r2Max-width/2/cos(endWedge);
                    end
                    borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rEnd/r2Max*cosEnd));
                    [xStart,yStart] = pol2cart(borderAngle,r2Max);
                    startWedge = borderAngle-lp/2*th2Moy+pi/2;
                end
                if rEnd >= r2Max
                    if rStart >=  r2Max-width/2/cos(startWedge)
                        width = (r2Max-rStart)*cos(startWedge)+width/2;
                        rStart = r2Max-width/2/cos(startWedge);
                    end
                    borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rStart/r2Max*cosStart));
                    [xEnd,yEnd] = pol2cart(borderAngle,r2Max);
                    endWedge = borderAngle-lp/2*th2Moy+pi/2;
                end
                if rStart <= r2Min
                    if rEnd <=  r2Min+width/2/cos(endWedge)
                        width = (rEnd-r2Min)*cos(endWedge)+width/2;
                        rEnd = r2Min+width/2/cos(endWedge);
                    end
                    borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rEnd/r2Min*cosEnd));
                    [xStart,yStart] = pol2cart(borderAngle,r2Min);
                    startWedge = borderAngle-lp/2*th2Moy+pi/2;
                end
                if rEnd <= r2Min
                    if rStart <= r2Min+width/2/cos(startWedge)
                        width = (rStart-r2Min)*cos(startWedge)+width/2;
                        rStart = r2Min+width/2/cos(startWedge);
                    end
                    borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rStart/r2Min*cosStart));
                    [xEnd,yEnd] = pol2cart(borderAngle,r2Min);
                    endWedge = borderAngle-lp/2*th2Moy+pi/2;
                end
                length = sqrt((xStart-xEnd)^2+(yStart-yEnd)^2);
                xPos = (xStart+xEnd)/2;
                yPos = (yStart+yEnd)/2;
                if drawLines == 1, PatchWedgedBlock
                    if nbQuadrants >= 2, PatchWedgedBlockQ2
                    end
                    if nbQuadrants == 4, PatchWedgedBlockQ4
                    end
                end
            end
        end
    end
    if rMin > rMax, outOfBounds = 1;
    end
end
toc
disp(sprintf('\n Number of discontinuities = %d.\n',jj))





% Graphics
la_tmp=.05;%
%line([-rMax rMax]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
%line([-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
%line([0 0],[-rMax rMax]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
%line([rMax/sqrt(2) -rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(rMax*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-rMax*(1+5*la_tmp),0,'$180^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),rMax/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),-rMax/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
text(-rMax*la_tmp,rMax*(1+2*la_tmp),'$90^{\circ}$')
text(-rMax*la_tmp,-rMax*(1+2*la_tmp),'$270^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),rMax/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),-rMax/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
xlabel('$\mu m$')
ylabel('$\mu m$')
set(gca,'xLim',[0 rMax]*(1+2*la_tmp))
set(gca,'ylim',[0 rMax]*(1+2*la_tmp))
if nbQuadrants >= 2, set(gca,'xLim',[-rMax rMax]*(1+2*la_tmp))
end
if nbQuadrants == 4, set(gca,'yLim',[-rMax rMax]*(1+2*la_tmp))
end
str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',periodStart,depth,periodRat,fillStart);
title(str,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')
%tick2latex

%save
print('-depsc2',sprintf('hybrid_periodRatio=%3.2f_nThLines=%d_offset=%d_rMax=%d.eps',periodRat,nThLines,offset,rMax), '-r300');


