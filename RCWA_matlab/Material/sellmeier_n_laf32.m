%Equation de Sellmeier n-Laf32

function[n]=sellmeier_n_laf32(l,K1,L1,K2,L2,K3,L3)

n=1+K1*l^2/(l^2-L1)+K2*l^2/(l^2-L2)+K3*l^2/(l^2-L3);