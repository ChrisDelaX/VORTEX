function[Eprof]=layerFab(F,E1,E2)

% profil blasé:
%
%   0=0%                                  1024=100%
%    |||||||||||                              |
%    |||||||||||||||||||||                    |
%    |||||||||||||||||||||||||||||||          |
%    ||||||||||||||||||||||||||||||||||||||||||     1ere couche

F1=F*1024;

for i=1:1024
    if i<=F1
        Eprof(1,i)=E2;
    else
        Eprof(1,i)=E1;
    end
end