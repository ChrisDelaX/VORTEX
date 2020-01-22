function[Eprof]=layer(F,E1,E2)

% profil rectangulaire simple:
%   0=0%                                  1024=100%
%    |             ||||||||||||||             |
%    |             F1||||||||||F2             |
%    |_____________||||||||||||||_____________|     1ere couche
%
% profil trapézoidal simple:
%   0=0%                                  1024=100%
%    |                ||||||                  |
%    |             F1||||||||F2               |
%    |__________||||||||||||||||||____________|     1ere couche


F1=(1-F)*1024/2;
F2=(1+F)*1024/2;

for i=1:1024
    if i<=F1
        Eprof(1,i)=E1;
    elseif i<=F2
        Eprof(1,i)=E2;
    else
        Eprof(1,i)=E1;
    end
end