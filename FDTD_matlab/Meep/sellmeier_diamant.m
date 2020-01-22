%Equation de Sellmeier diamant

function[n]=sellmeier_diamant(l,A,B,C,D,E)

l=l*1e3;

n=(A+(B*(l^2)/((l^2)-C^2))+(D*(l^2)/((l^2)-E^2)));