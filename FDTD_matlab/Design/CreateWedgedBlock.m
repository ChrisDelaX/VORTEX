function[x,y]=CreateWedgedBlock(xPos,yPos,width,length,startWedge,endWedge,twist,blockColor,nbQuadrants)

% % TEST
% close all
% newFig
% axis equal
% xPos = 0;
% yPos = 0;
% width = 7;
% length = 30;
% twist = deg2rad(10);
% startWedge = deg2rad(20);
% endWedge = deg2rad(-45);

lengthTempStart = width/2*tan(startWedge);
lengthTempEnd = width/2*tan(endWedge);
[RotateX,RotateY] = RotateAngle(width/2,length/2+lengthTempEnd,twist);
x(1) = xPos + RotateX;
y(1) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(width/2,-length/2+lengthTempStart,twist);
x(2) = xPos + RotateX;
y(2) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(-width/2,-length/2-lengthTempStart,twist);
x(3) = xPos + RotateX;
y(3) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(-width/2,length/2-lengthTempEnd,twist);
x(4) = xPos + RotateX;
y(4) = yPos + RotateY;
patch(x,y,blockColor,'edgecolor','none')
if nbQuadrants >= 2, patch(-x,y,blockColor,'edgecolor','none')
end
if nbQuadrants ==4    
patch(-x,-y,blockColor,'edgecolor','none')
patch(x,-y,blockColor,'edgecolor','none')
end
