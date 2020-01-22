%Equation de Sellmeier Si

function[n]=sellmeier_znseb(l,A,B,C,D,E)

n=(A+B*l+C*l^2+D*l^3)^2;