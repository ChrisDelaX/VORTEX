% by Christian Delacroix (2014)

fprintf(fid,meepMakeBlock(nmat,supCenter,supSize,e1,e2,e3));

%sub createHWP(cx,cy,cz,vx,vz,rotz,swa,nb45,xSize,ySize,period)
    dim diag,orient1,orient2,xMid,yMid,r,r2,rPos,dxLength,dyLength,dWedge,inp,np,xPos,yPos
    twist = modulo(rotz,pi)
    diag = atn(ySize/xSize)
    orient1 = pi/2-twist
    orient2 = abs(orient1)
    if orient2 >= diag then orient2 = pi/2-orient2
    xMid = xSize/2-ySize/2*tan(orient2)
    yMid = ySize/2-xSize/2*tan(orient2)    
    dxLength = ySize/cos(orient2)
    dyLength = xSize/cos(orient2)
    r = period/cos(orient2)   
    if abs(orient2) > 0 then r2 = period/sin(2*orient2)
    dWedge = orient2*sign(orient1)
    if abs(orient1) >= diag then '' ALONG X-AXIS
        if abs(orient1) = diag then 'central line
            call createTrapezoidWedged(cx,cy,cz,vx,dxLength,vz,twist,swa,nb45,0,0)
        else
            call createTrapezoidWedged(cx,cy,cz,vx,dxLength,vz,twist,swa,nb45,-dWedge,-dWedge) 
        end if
        np = floor((xSize*cos(orient2)+ySize*sin(orient2)-period+vx)/2/period)
        for inp=1 to np
            xPos = r*inp
            if xPos <=
 xMid then 'longer lines
                call createTrapezoidWedged(cx+xPos,cy,cz,vx,dxLength,vz,twist,swa,nb45,-dWedge,-dWedge)
                call createTrapezoidWedged(cx-xPos,cy,cz,vx,dxLength,vz,twist,swa,nb45,-dWedge,-dWedge)
            else 'shorter lines
                if abs(orient2) > 0 then dxLength = (xSize-xMid-r*inp)/sin(orient2)
                rPos = r2*(inp-xMid/r)
                xPos = rPos*sin(orient2)+xMid
                yPos = rPos*cos(orient2)*sign(orient1)
                dWedge = -orient2*sign(orient1) 
                call createTrapezoidWedged(cx+xPos,cy+yPos,cz,vx,dxLength,vz,twist,swa,nb45,dWedge,dWedge+pi/2)
                call createTrapezoidWedged(cx-xPos,cy-yPos,cz,vx,dxLength,vz,twist,swa,nb45,dWedge+pi/2,dWedge)
            end if
        next
    else '' ALONG Y-AXIS
        call createTrapezoidWedged(cx,cy,cz,vx,dyLength,vz,twist,swa,nb45,dWedge,dWedge) 'central line
        np = floor((ySize*cos(orient2)+xSize*sin(orient2)-period+vx)/2/period)
        for inp=1 to np
            yPos = r*inp
            if yPos <= yMid then 'longer lines
                call createTrapezoidWedged(cx,cy+yPos,cz,vx,dyLength,vz,twist,swa,nb45,dWedge,dWedge)
                call createTrapezoidWedged(cx,cy-yPos,cz,vx,dyLength,vz,twist,swa,nb45,dWedge,dWedge)
            else 'shorter lines
                if abs(orient2) > 0 then dyLength = (ySize-ymid-r*inp)/sin(orient2)
                rPos = r2*(inp-yMid/r)
                xPos = rPos*cos(orient2)*sign(orient1)
                yPos = rPos*sin(orient2)+yMid
                dWedge = orient2*sign(orient1) + pi/2*(1-sign(orient1))/2
                call createTrapezoidWedged(cx-xPos,cy-yPos,cz,vx,dyLength,vz,twist,swa,nb45,dWedge,dWedge+pi/2)
                call createTrapezoidWedged(cx+xPos,cy+yPos,cz,vx,dyLength,vz,twist,swa,nb45,dWedge+pi/2,dWedge)
            end if
        next
    end if
end sub