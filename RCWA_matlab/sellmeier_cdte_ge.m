%Equation de Sellmeier CdTe et Ge

function[n]=sellmeier_cdte_ge(l,A,B,C,D,E)

n=(A+(B*(l^2)/((l^2)-C))+(D*(l^2)/((l^2)-E)));