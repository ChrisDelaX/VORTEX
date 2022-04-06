function[Eprof]=layer9(F,Etmp,E1,E2)

% profil rectangulaire simple avec couche d'arret discontinue
% remplissage ?????????? période indéfinie !
%
%   0=0%                                       1024=100%
%    ||||||F1     F2||||||||||||F3     F4|||||||||||
%    |aaaaaF1aaaaaF2||||||||||||F3aaaaaF4aaaaaaaaaa|

F1=1/2*F*1024/2;
F2=2/2*F*1024/2;
F3=4/2*F*1024/2;
F4=5/2*F*1024/2;

for i=1:1024
    if i<=F1
        Eprof(1,i)=Etmp;
    elseif i<=F2
        Eprof(1,i)=E1;
    elseif i<=F3
        Eprof(1,i)=E2;
    elseif i<=F4
        Eprof(1,i)=E1;
    else
        Eprof(1,i)=Etmp;
    end
end