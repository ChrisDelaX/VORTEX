function disk(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
mygrey = [0.35,0.35,0.35];
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
patch(x+xp,y+yp,mygrey,'edgecolor',mygrey)
end