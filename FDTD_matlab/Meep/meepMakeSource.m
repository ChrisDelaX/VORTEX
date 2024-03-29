% by Christian Delacroix (2014)

function str = meepMakeSource(dx,dy,dfreq,Ename,Einput)

str = ['(make source (src (make continuous-src ' ...
    '(frequency ' num2str(dfreq) ')))' ...
    '(component ' Ename ')' ...
    '(center 0 0 ' num2str(Einput) ')' ...
    '(size ' num2str(dx) ' ' num2str(dy) ' 0))\n'];
