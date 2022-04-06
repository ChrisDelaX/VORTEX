function[Eprof]=layer12(F,Etmp,E1,E2)

% profil rectangulaire avec couche d'arret discontinue et antireflet continue
% remplissage ?????????? période indéfinie !
%
%   0=0%                                                    1024=100%
%    |RRRRRF1RRRRRRRRRRRRRRF2RRRRRRRRRRRRF3RRRRRRRRRRRRRRF4RRRRR|
%    ||||||F1              F2||||||||||||F3              F4||||||
%    |aaaaaF1aaaaaaaaaaaaaaF2||||||||||||F3aaaaaaaaaaaaaaF4aaaaa|

F1=(F/2)*1024/2;
F2=(1-F/2)*1024/2;
F3=(1+F/2)*1024/2;
F4=(2-F/2)*1024/2;

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