function[n]=sellmeier_lah83(l,A1,B1,A2,B2,A3,B3)
% l=l*1000;

n=1+A1*l^2/(l^2-B1)+A2*l^2/(l^2-B2)+A3*l^2/(l^2-B3);