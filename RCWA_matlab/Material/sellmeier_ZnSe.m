%Equation de Sellmeier ZnSe

function[n]=sellmeier_ZnSe(l,A,B,C,D,E,F,G)

n=(A+(B*(l^2)/((l^2)-C))+(D*(l^2)/((l^2)-E))+(F*(l^2)/((l^2)-G)));