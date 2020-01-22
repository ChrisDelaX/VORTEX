function[RotateX,RotateY] = RotateAngle(dX,dY,dAngle)
    RotateX = dX*cos(dAngle) - dY*sin(dAngle);
    RotateY = dX*sin(dAngle) + dY*cos(dAngle);