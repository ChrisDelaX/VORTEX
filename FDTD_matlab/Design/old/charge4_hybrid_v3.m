clear all
close all

% INPUT PARAMETERS
% ----------------
lp = 4;     % topological charge
periodStart = 1.42;%1;%
periodEnd = 0.95*periodStart;
fillStart = 0.46; %0.2;%
fillEnd = 0.49; %0.5;%
nPeriods = 20;%60;%120;%40;%
nRLines = 1;%2;% % lines radius discretization (number)
nThLines = 4;%16;%  % lines theta discretization (number)
offset = 3;%2;%  % curves theta discretization (in periods!!!)

% Graphics, numerics
% ------------------
newFig
curveCellColor = myblue;
lineCellColor = mygreen;
blockColor = mygrey;%'k';%
xCenter = 0;
yCenter = 0;
maxIter = 1e2;
thDisc = deg2rad(3);%5 % curves precision (in rad)
nbQuadrants = 4;  % 1,2 or 4
drawCells = 0;  % 0 or 1
drawCurves = 0;  % 0 or 1
drawLines = 1;  % 0 or 1


% Special cases
% -------------
%offset = 0;                 % only curves!
periodEnd = periodStart;    % only lines!
fourCases = 4;%1;% 
for ff = 1:fourCases
if fourCases == 4
nbQuadrants = 1;
%nThLines = ff*4+4;%2;%
nRLines = ff*2;
end
    

% Merging verification
% --------------------
depth = 5.3;% %grating depth
swa = deg2rad(0);   % sidewall angle
if periodEnd < 2*depth*tan(swa)/(1-fillEnd)
    2*depth*tan(swa)/(1-fillEnd)
    error('periodRat too small. Merging of the walls.')
end


% Initializing figure and values
% ------------------------------
widthStart = fillStart*periodStart;
widthEnd = fillEnd*periodEnd;
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
nThLinestmp = nThLines;
rMax = nPeriods*periodStart;
rMin = 0;
iter=0;
tic
datestr(now)
while iter < maxIter   % theta discretization for curves
    
    % Curves along x-axis
    % -------------------
    nbStart = ceil(rMin/stepThStart+(1-fillStart)/2);
    nbEnd = ceil(rMax/stepThEnd-(1-fillStart)/2);
    for iR=nbStart:nbEnd
        % Draws Cells
        if iR == nbStart
            [xCell1(1),yCell1(1)] = pol2cart(thStart,rMax);
            [xCell1(2),yCell1(2)] = pol2cart(thStart,rMin);
            [xCell1(3),yCell1(3)] = pol2cart(thEnd,rMin);
            [xCell1(4),yCell1(4)] = pol2cart(thEnd,rMax);
        end
        % Draws Curves
        rPhi = stepPhi*(2*iR-1)/4;
        rThStart = 2*rPhi*cos(thStart);
        rThEnd = 2*rPhi*cos(thEnd);
        edgeStart = 0;
        edgeEnd = 0;
        thStart0 = thStart;
        thEnd0 = thEnd;
        widthStart0 = widthStart;
        widthEnd0 = widthEnd;
        if rThStart > rMax
            edgeStart = 1;
            thStart0 = acos(rMax/(2*rPhi));
            widthStart0 = widthStart+(widthEnd-widthStart)*(thStart0-thStart)/(thEnd-thStart);
        end
        if rThEnd < rMin
            edgeEnd = 1;
            thEnd0 = acos(rMin/(2*rPhi));
            widthEnd0 = widthEnd+(widthStart-widthEnd)*(thEnd0-thEnd)/(thStart-thEnd);
        end
        phiStart0 = 2*thStart0;
        phiEnd0 = 2*thEnd0;
        startWedge0 = thStart0+edgeStart*pi/2;
        endWedge0 = thEnd0+edgeEnd*pi/2;
        nThCurves = ceil((phiEnd0-phiStart0)/thDisc);
        for iTH=1:nThCurves
            width = widthStart0+(widthEnd0-widthStart0)/nThCurves*(iTH-1/2);
            phiStart = phiStart0+(phiEnd0-phiStart0)/nThCurves*(iTH-1);
            phiEnd = phiStart0+(phiEnd0-phiStart0)/nThCurves*(iTH);
            length = 2*rPhi*sin((phiEnd-phiStart)/2);
            twist = (phiStart+phiEnd)/2;
            [xPos,yPos] = pol2cart(twist,rPhi);
            xPos = xCenter + xPos + rPhi;
            yPos = yCenter + yPos;
            startWedge = startWedge0+(endWedge0-startWedge0)/nThCurves*(iTH-1)-twist;
            endWedge = startWedge0+(endWedge0-startWedge0)/nThCurves*(iTH)-twist;
            if drawCurves == 1
                if ff == 1, CreateWedgedBlock(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)                  
                elseif ff == 2, CreateWedgedBlockQ2(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                elseif ff == 3, CreateWedgedBlockQ3(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)  
                elseif ff == 4, CreateWedgedBlockQ4(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                end
            end
        end
    end
    % Increments
    iter=iter+1;
    cosStart = periodRat^(2*iter/lp);
    cosEnd = periodRat^(2*(iter+1)/lp);
    thStart = acos(cosStart);
    thEnd = acos(cosEnd);
    thMoy = (thStart+thEnd)/2;
    stepPhi = periodStart/periodRat^(iter);
    stepThStart = stepPhi*cos(thStart);
    stepThEnd = stepPhi*cos(thEnd);
    r2Min = rMin;
    rMin = offset/2/tan(thMoy-thStart);
    r2Max = min(rMin,rMax);
    
    
    % Lines along y-axis
    % ------------------
    nl = 1;
    rtmp = (r2Max-r2Min)/nRLines;
    r2Maxtmp = r2Max;
    r2Max = r2Min+rtmp;
    nThLines = nThLinestmp;
    while nl <= nRLines
        theta = (pi/2-thStart)/nThLines;
        for iTH=1:nThLines
            % Draws Cells
            th2Start = thStart+(pi/2-thStart)/nThLines*(iTH-1);
            th2End = min(pi/2,thStart+(pi/2-thStart)/nThLines*(iTH));
            th2Moy = (th2Start+th2End)/2;
            cosEnd = cos(th2End);
            cosStart = cos(th2Start);
            [xCell2(1),yCell2(1)] = pol2cart(th2End,r2Max);
            [xCell2(2),yCell2(2)] = pol2cart(th2End,r2Min);
            [xCell2(3),yCell2(3)] = pol2cart(th2Start,r2Min);
            [xCell2(4),yCell2(4)] = pol2cart(th2Start,r2Max);
            if drawCells == 1, 
                createCell(xCell2,yCell2,lineCellColor,nbQuadrants);
                createCell(xCell1,yCell1,curveCellColor,nbQuadrants);
            end
            % Draws Lines
            twist = lp/2*th2Moy;
            step2 = periodStart/cos(th2Moy);
            step2End = periodStart/cos(th2End);
            nbStart2 = floor(r2Min/(cos(th2Moy-th2Start)+sin(th2Moy-th2Start)*tan(th2End)) /step2 +(1-fillStart)/2);
            nbEnd2 = floor(r2Max*(cos(th2Moy-th2Start)+sin(th2Moy-th2Start)*tan(th2Moy)) /step2 -(1-fillStart)/2);
            iR=0;
            for iR=nbStart2:nbEnd2
                width = fillStart*periodStart;
                r = (iR+1/2)*step2 ;
                rStart = r*cos(th2Moy)/cosEnd;
                rEnd = r*cos(th2Moy)/cosStart;
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
                        borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rEnd/r2Max*cosStart));
                        [xStart,yStart] = pol2cart(borderAngle,r2Max);
                        startWedge = borderAngle-lp/2*th2Moy+pi/2;
                    end
                    if rEnd >= r2Max
                        if rStart >=  r2Max-width/2/cos(startWedge)
                            width = (r2Max-rStart)*cos(startWedge)+width/2;
                            rStart = r2Max-width/2/cos(startWedge);
                        end
                        borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rStart/r2Max*cosEnd));
                        [xEnd,yEnd] = pol2cart(borderAngle,r2Max);
                        endWedge = borderAngle-lp/2*th2Moy+pi/2;
                    end
                    if rStart <= r2Min
                        if rEnd <=  r2Min+width/2/cos(endWedge)
                            width = (rEnd-r2Min)*cos(endWedge)+width/2;
                            rEnd = r2Min+width/2/cos(endWedge);
                        end
                        borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rEnd/r2Min*cosStart));
                        [xStart,yStart] = pol2cart(borderAngle,r2Min);
                        startWedge = borderAngle-lp/2*th2Moy+pi/2;
                    end
                    if rEnd <= r2Min
                        if rStart <= r2Min+width/2/cos(startWedge)
                            width = (rStart-r2Min)*cos(startWedge)+width/2;
                            rStart = r2Min+width/2/cos(startWedge);
                        end
                        borderAngle = pi+(lp/2)*th2Moy+(pi/2+asin(rStart/r2Min*cosEnd));
                        [xEnd,yEnd] = pol2cart(borderAngle,r2Min);
                        endWedge = borderAngle-lp/2*th2Moy+pi/2;
                    end
                    length = sqrt((xStart-xEnd)^2+(yStart-yEnd)^2);
                    xPos = xCenter + (xStart+xEnd)/2;
                    yPos = yCenter + (yStart+yEnd)/2;
                    if drawLines == 1
                        if ff == 1, CreateWedgedBlock(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                        elseif ff == 2, CreateWedgedBlockQ2(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                        elseif ff == 3, CreateWedgedBlockQ3(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                        elseif ff == 4, CreateWedgedBlockQ4(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)
                        end
                    end
                end
            end
        end
        % Increments
        nl = nl+1;
        r2Min = r2Min+rtmp;
        r2Max = r2Max+rtmp;
        nThLines = nThLines+nThLinestmp;
    end
    nThLines = nThLinestmp;
    if rMin > rMax, break
    end
end
end
toc
disp(sprintf('\n Number of curve discontinuities = %d.\n',iter))



% Graphics
% --------
grid off
set(gca,'visible','off')
axis equal
xlabel('$\mu m$')
ylabel('$\mu m$')
la_tmp=.05;%
set(gca,'xLim',[0 rMax]*(1+2*la_tmp))
set(gca,'ylim',[0 rMax]*(1+2*la_tmp))
if nbQuadrants >= 2 || fourCases == 4, set(gca,'xLim',[-rMax rMax]*(1+2*la_tmp))
end
if nbQuadrants == 4 || fourCases == 4, set(gca,'yLim',[-rMax rMax]*(1+2*la_tmp))
end
line([-rMax rMax]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
line([-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
line([0 0],[-rMax rMax]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
line([rMax/sqrt(2) -rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(rMax*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-rMax*(1+5*la_tmp),0,'$180^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),rMax/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),-rMax/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
text(-rMax*la_tmp,rMax*(1+2*la_tmp),'$90^{\circ}$')
text(-rMax*la_tmp,-rMax*(1+2*la_tmp),'$270^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),rMax/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),-rMax/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',periodStart,depth,periodRat,fillStart);
title(str,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')
%tick2latex

%save
%print('-depsc2',sprintf('01_lines_rMax=%2.0f_nThLines=%d.eps',rMax,nThLinestmp), '-r300');
%print('-depsc2',sprintf('02_curves_rMax=%2.0f_periodRatio=%3.2f.eps',rMax,periodRat), '-r300');
print('-depsc2',sprintf('03_rings_rMax=%2.0f_nThLines=%d_nRLines=%d.eps',rMax,nThLinestmp,nRLines), '-r300');
%print('-depsc2',sprintf('04_hybrid-cells_rMax=%2.0f_periodRatio=%3.2f_nThLines=%d.eps',rMax,periodRat,nThLinestmp), '-r300');
%print('-depsc2',sprintf('05_hybrid_rMax=%2.0f_periodRatio=%3.2f_nThLines=%d.eps',rMax,periodRat,nThLinestmp), '-r300');


