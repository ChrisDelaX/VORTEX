% by Christian Delacroix (2014)

fid = fopen([sgvc '.ctl'],'w'); % Meep Control Script (.ctl)
%fprintf(fid,'(set! print-ok? false)\n');

%% Settings, pml layers, sources
% ******************************
fprintf(fid,['(set! resolution ' num2str(rez) ')\n']);
if force_complex fprintf(fid,['(set! force-complex-fields? true)\n']); end
fprintf(fid,['(set! pml-layers (list\n' ['   ' meepMakePML(1)] '))\n']);
fprintf(fid,'(set! sources (list\n');
switch pola
    case 'LHC'
        fprintf(fid,['   ' meepMakeSource(xArea,yArea,1/wavelength,'Ex',zInput)]);
        fprintf(fid,['   ' meepMakeSource(xArea,yArea,1/wavelength,'Ey',zInput-wavelength/4)]);
    case 'RHC'
        fprintf(fid,['   ' meepMakeSource(xArea,yArea,1/wavelength,'Ex',zInput-wavelength/4)]);
        fprintf(fid,['   ' meepMakeSource(xArea,yArea,1/wavelength,'Ey',zInput)]);
end
fprintf(fid,'))\n');

%% Set geometry
% *************
e1 = [cosd(yTilt) 0 -sind(yTilt)];
e2 = [0 1 0];
e3 = [sind(yTilt) 0 cosd(yTilt)];
fprintf(fid,['(set! geometry-lattice (make lattice ' ...
    '(size ' num2str(xArea) ' ' num2str(yArea) ' ' num2str(zArea) ')))\n']);
fprintf(fid,'(set! geometry (list\n');
switch sgvc
    case 'AGPM'
        for i=1:nPeriods
            rTop = width/2+(nPeriods-i+1)*period;
            fprintf(fid,['   ' meepMakeCone(nmat,center,e3,depth,rTop+TBdiff,rTop)]);
            rTop = rTop-width;
            fprintf(fid,['   ' meepMakeCone(1,center,e3,depth,rTop-TBdiff,rTop)]);
        end
        fprintf(fid,['   ' meepMakeCone(nmat,center,e3,depth,width/2+TBdiff,width/2)]);
    case 'multi'
        xSize = xArea/5;%2
        ySize = yArea;%/5'2
        dp = 0.02;
        for i=-2:1:2
            for j=1:ceil(yArea/period)
                hwpC = [xOffset+i*xSize yOffset-yArea/2+(j-1/2)*period zCenter];
                hwpS = [xSize width+i*dp depth];
                fprintf(fid,['   ' meepMakeBlock(nmat,hwpC,hwpS,e1,e2,e3)]);
            end
        end
end
supCenter = [zCenter-sind(yTilt)*(depth/2+zSupport/2) 0 zCenter-cosd(yTilt)*(depth/2+zSupport/2)];
supSize = [2*xArea,2*yArea,zSupport];
fprintf(fid,['   ' meepMakeBlock(nmat,supCenter,supSize,e1,e2,e3)]);
fprintf(fid,'))\n\n');


%% set boundaries and output fields
% *********************************
fprintf(fid,['(run-until ' num2str(steps) '\n']);
fprintf(fid,['   (at-beginning\n' ...
    '      output-epsilon\n' ...
    '   )\n']);
    %'(at-end output-efield-x)\n'...
    %'(at-end output-efield-y)\n'...
    

for i = zOutput %for writing out fields in specified planes
    fprintf(fid,['   (at-end\n' ...
        '      (in-volume \n' ...
        '         (volume (center 0 0 ' num2str(zCenter-depth/2-i) ') (size ' num2str(xArea*(1-crop)) ' ' num2str(yArea*(1-crop)) ' 0))\n' ...
        '         (to-appended "Ex_z=' num2str(i) '" output-efield-x)\n' ...
        '         (to-appended "Ey_z=' num2str(i) '" output-efield-y)\n' ...
        '      )\n' ...
        '   )\n']);
end

%% Time Series
% ************
if rec_time_series %for writing out a field in a plane at every timestep after recT
    fprintf(fid,['(after-time ' num2str(recT) ' (in-volume \n(volume (center 0 0 ' num2str(zCenter-depth/2-recZ) ') '...
        '(size ' num2str(xArea*(1-crop)) ' ' num2str(yArea*(1-crop)) ' 0))\n (to-appended "Ex-xyt" output-efield-x)(to-appended "Ey-xyt" output-efield-y)))']);
end

%% end
% ****
fprintf(fid,')\n');
fclose(fid);