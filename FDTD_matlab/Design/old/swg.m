clear all
close all

lp = 4;
period = 1.4;%'1.42'
depth = 5.2;%'4.7'
fill = 0.45;%'0.41'
%swAngle = 3;%'6'
width = period*fill;%
%nWidths = 10;%'5'&
nThetas = 16;%24;%64;%   ''multiple of 4
nPeriods = 10;%30;%7;%
orient = 1*pi/6;%-4;%
rMax = period*(nPeriods+0.5);

%center
x0 = 10;
y0 = 15;

% Figures
% -------
newFig
grid off
%set(gca,'visible','off')
axis equal
winSize = 1.1;
axis([-rMax*winSize+x0 rMax*winSize+x0 -rMax*winSize+y0 rMax*winSize+y0])


orient = mod(orient,pi);
orient2 = mod(orient,pi/2);
orient3 = pi/4-abs(orient2-pi/4);
orient4 = pi/2-orient3;
twist = -orient;
lgth = 2*rMax/cos(orient3);
rmid = rMax*(1-1/tan(orient4));
r = period/cos(orient3);
r2 = period/cos(orient4-orient3);
nPeriods2 = ceil((rMax/period-0.5)*(1/cos(orient3)+(1-tan(orient3))*sin(orient3)));


% central line
% ------------
length = lgth;
if abs(orient-pi/2) >= pi/4
    startWedge = -orient3*sign(orient-pi/2);
    endWedge = -orient3*sign(orient-pi/2);
else
    startWedge = orient3*sign(orient-pi/2);
    endWedge = orient3*sign(orient-pi/2);
end
xPos = x0;
yPos = y0;
twist = -orient;
CreateTrapezoidWedged

tic
depart_tmp=now;
depart=datestr(depart_tmp)
count = 0;
for ii=1:nPeriods2    
    rPos = r*ii;
    % longer lines
    % ------------
    if rPos <= rMax*(1-tan(orient3))
        count = count + 1;
        length = lgth;
        % along x-axis
        if abs(orient-pi/2) >= pi/4
            startWedge = -orient3*sign(orient-pi/2);
            endWedge = -orient3*sign(orient-pi/2);
            xPos = x0+rPos;
            yPos = y0;
            CreateTrapezoidWedged
            xPos = x0-rPos;
            yPos = y0;
            CreateTrapezoidWedged
        % along y-axis
        else 
            startWedge = orient3*sign(orient-pi/2);
            endWedge = orient3*sign(orient-pi/2);
            yPos = y0+rPos;
            xPos = x0;
            CreateTrapezoidWedged
            yPos = y0-rPos;
            xPos = x0;
            CreateTrapezoidWedged
        end
    % shorter lines
    % -------------
    else
        r2Pos = r2*(ii-count-(rmid-(r*count))/r);
        length = lgth/2 - (rPos-rMax)/sin(orient3);
        % along x-axis
        if abs(orient-pi/2) >= pi/4
            startWedge = -pi/2-orient3*sign(orient-pi/2);
            endWedge = -orient3*sign(orient-pi/2);
            xPos = x0+rmid+r2Pos*cos(orient4);
            yPos = y0+r2Pos*sin(orient4)*sign(orient-pi/2);
            CreateTrapezoidWedged
            endWedge = -pi/2-orient3*sign(orient-pi/2);
            startWedge = -orient3*sign(orient-pi/2);
            xPos = x0-rmid-r2Pos*cos(orient4);
            yPos = y0-r2Pos*sin(orient4)*sign(orient-pi/2);
            CreateTrapezoidWedged
        % along y-axis
        else
            startWedge = orient3*sign(orient-pi/2) - pi/2*(1+sign(orient-pi/2))/2;
            endWedge = orient3*sign(orient-pi/2) - pi/2*(1-sign(orient-pi/2))/2;
            yPos = y0-rmid-r2Pos*cos(orient4);
            xPos = x0-r2Pos*sin(orient4)*sign(orient-pi/2);
            CreateTrapezoidWedged
            endWedge = orient3*sign(orient-pi/2) - pi/2*(1+sign(orient-pi/2))/2;
            startWedge = orient3*sign(orient-pi/2) - pi/2*(1-sign(orient-pi/2))/2;
            yPos = y0+rmid+r2Pos*cos(orient4);
            xPos = x0+r2Pos*sin(orient4)*sign(orient-pi/2);
            CreateTrapezoidWedged
        end
    end    
end


%save
print('-depsc2',sprintf('swg_rMax=%3.2f_orient=%3.2f.eps',rMax,orient), '-r300');









