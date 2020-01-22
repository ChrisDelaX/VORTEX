% % TEST
% close all
% newFig
% axis equal
% xPos = 0;
% yPos = 0;
% width = 7;
% period = 10;
% length = 30;
% twist = deg2rad(0);
% startWedge = deg2rad(50);
% endWedge = deg2rad(-45);

lengthTempStart = width/2*tan(startWedge);
lengthTempEnd = width/2*tan(endWedge);
lengthTemp = length + (lengthTempStart-lengthTempEnd);
[RotateX,RotateY] = RotateAngle(width/2,length/2+lengthTempStart,twist);
x(1) = xPos + RotateX;
y(1) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(width/2,-length/2+lengthTempEnd,twist);
x(2) = xPos + RotateX;
y(2) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(-width/2,-length/2-lengthTempEnd,twist);
x(3) = xPos + RotateX;
y(3) = yPos + RotateY;
[RotateX,RotateY] = RotateAngle(-width/2,length/2-lengthTempStart,twist);
x(4) = xPos + RotateX;
y(4) = yPos + RotateY;
patch(x,y,mygrey)%,'edgealpha',0)
%patch(x,y,[1-1/nLevels*nl,0.2,0+1/nLevels*nl])%,'edgealpha',0)






