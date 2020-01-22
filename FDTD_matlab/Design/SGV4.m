clear all
close all

% Output files
% ------------
filename = 'test_charge4';
fid = fopen([filename '.scr'],'w');
prec = 8; % floating point precision

% INPUT PARAMETERS
% ----------------
nPeriods = 100;
lp = 4;             % topological charge
periodStart = 1.42; %
periodEnd = 0.8*periodStart;
fillStart = 0.588;   % rcwa result (trapez=0.46)
fillEnd = 0.626;     % rcwa result (trapez=0.49)
depth = 4.558;      % grating depth
swa = deg2rad(0);   % sidewall angle
xCenter = 0;
yCenter = 0;
drawCurves = 1;  % 0 for NO curved lines, otherwise 1
drawLines = 1;  % 0 for NO straight lines, otherwise 1
drawCells = 0;  % 0 or 1


% DISCRETIZATION
% --------------
nRLines = 1;    % lines radius discretization (number)
nThLines = 6;%16;   % lines theta discretization (number)  = 2 for the EOPM, 6 for hybrid (with curves)
offset = 8;     % curves theta discretization (in periods!!!)
maxIter = 1e2;
thDisc = deg2rad(3);%5  % curves precision (in rad)
nbQuadrants = 4;%1;%  % 1,2 or 4


% Special case: 4 different quadrants
% -----------------------------------
fourCases = 1;%4;   % Leave 1, unless you specify different designs for each quadrant (then 4)
for ff = 1:fourCases
if fourCases == 4
    nbQuadrants = 1;
    nThLines = ff*4+4;%2;%
    %nRLines = ff*2;
end


% Initializing figure and values
% ------------------------------
% Merging verification
if periodEnd < 2*depth*tan(swa)/(1-fillEnd)
    2*depth*tan(swa)/(1-fillEnd)
    error('Merging of the walls: reduce sidewall angle.')
end
newFig
curveCellColor = myblue;
lineCellColor = mygreen;
blockColor = 'k';%mygrey;%
offset = offset * drawLines;                                        % only curved lines!
periodEnd = periodStart + (periodEnd-periodStart) * drawCurves;     % only straight lines!
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
                if ff == 1, 
                    [x,y] = CreateWedgedBlock(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants);
                    fprintf(fid,['pline ' num2str(x(1),prec) ',' num2str(y(1),prec) ' ' num2str(x(2),prec) ',' num2str(y(2),prec) ' ' num2str(x(3),prec) ',' num2str(y(3),prec) ' ' num2str(x(4),prec) ',' num2str(y(4),prec) ' c\n']);
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
                CreateCell(xCell2,yCell2,lineCellColor,nbQuadrants);
                CreateCell(xCell1,yCell1,curveCellColor,nbQuadrants);
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
                        if ff == 1, 
                            [x,y]=CreateWedgedBlock(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants);
                            fprintf(fid,['pline ' num2str(x(1),prec) ',' num2str(y(1),prec) ' ' num2str(x(2),prec) ',' num2str(y(2),prec) ' ' num2str(x(3),prec) ',' num2str(y(3),prec) ' ' num2str(x(4),prec) ',' num2str(y(4),prec) ' c\n']);
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
xlabel('$\mu m$')
ylabel('$\mu m$')
zoom = 150;%60;%10;%40;%
myaxis = [-zoom zoom -zoom zoom];
axis equal; axis(myaxis)


%save
print('-depsc2',[filename '.eps'], '-r300');


