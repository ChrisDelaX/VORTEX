%% The Subwavelength Grating Vortex Coronagraph (SGVC) simulator
% by Christian Delacroix (2014)
% using a FDTD freeware called MEEP
% http://ab-initio.mit.edu/wiki/index.php/Meep
% Dug out in Jan 2020, by Christian Delacroix & Lorenzo König

clear all, close all

%% Case selection, grating parameters
% ***********************************
sgvc = 'multi';                  % case select: HWP,AGPM,multi,Hybrid
lp = 2;%4;%10;%                    % topological charge
offsetCurves = 3;%2;%           % curves radius offset (in periods!!!)
nRLines = 1;%2;%                % number of line rings
nThLines = 4;%16;%              % number of line slices in 90°
nCurve = 15;%30;%               % number of blocks for a 90° curve
nWedge = 7;%10;%                % number of blocks for a 45° wedge
fourCases = 1;%4;%              % 1 or 4
nQuadrants = 4;%1;%             % 1,2 or 4
maxIter = 100;
swa = 0;%3;%                    % sidewall angle (deg)
depth = 5.15;%5;%4.8;%          % grating depth (µm)
TBdiff = tand(swa)*depth;       % width diff between top and bottom, with tilted sidewall
periodStart = 1.4;%1.42;%       % (µm)
periodEnd = 0.8*periodStart;    % (µm)
fillStart = 0.87/periodStart;%0.61;%0.46;%
fillEnd = fillStart-0.1*(periodEnd-periodStart); % empiric
period = periodStart;
fill = fillStart;
width = period*fill;
if periodEnd < 2*depth*tand(swa)/(1-fillEnd)	% merging verification
    error('PeriodRat too small. Merging of the walls.')
end
drawLines = 1;                  % 0 or 1
%offsetCurves = 0;                      % only curves!
drawCurves = 1;                 % 0 or 1
%periodEnd = (1-1e-6)*periodStart;      % only lines!
%quad =
%horiz =
%vertic =
%ringPixNum =


%% Wafer parameters
% *****************
rez = 8;%5;%2;%20;%             % resolution (points per µm)
xArea = 10;%                  % (µm)
yArea = 5;                  % (µm)
zArea = 5.25;%20;                  % (µm)
xnpts = ceil(xArea*rez);                     % even number of points along x-axis
ynpts = ceil(yArea*rez);                          % even number of points along y-axis
znpts = ceil(zArea*rez);                        % even number of points along z-axis
nPeriods = ceil((xArea*sqrt(2)/2-period/2)/period);
xOffset = 0;%-25.6/2;%                  % grating center (µm)
yOffset = 0;%-25.6/2;%                  % grating center (µm)
zOffset = 0;%6.85;                      % air before the grating (min 2µm)
yTilt = 0;                              % tilt around y axis (deg)
zMargin = 20;%                          % substrate after the grating (no limit)
zSupport = zArea-zOffset-depth+zMargin;
zCenter = zArea/2-zOffset-depth/2;
center = [xOffset yOffset zCenter];     % grating center

%% Input wave
% ***********
zFields = false;                      % 'true' not yet implemented
force_complex = true;%false;%       % use complex fields
steps = 200;%300;                            % number of time steps
amplitude = 1;
wavelength = 3.8;                       % (µm)
nmat = sqrt(sellmeier_diamant(wavelength,1,0.3306,175.0,4.3356,106.0)); % refractive index (2.38)
pola = 'RHC';                           % RHC, LHC, VL, HL
zInput = zArea/2-1;%0.245;%             % input plane z position (µm), value from OPTIWAVE = 0.245


%% Observation planes
% *******************
% save E-fields at specific planes AFTER the grating (µm)
%zOutput = [-depth -depth/2 0 zSupport-zMargin];%                
zOutput = [0];% zSupport-zMargin];%                
zOutputAll = false;                     % save complete Exy(z)
crop = 0;%0.05;%                        % crop margin (fraction)
rec_time_series = false; %output a time series of a field component in a Z-plane.
recT = 290; %from this time step, start recording time series data.
recZ = 13; %how far under the AGPM to record data (plane not tilted with AGPM).

%% Run Meep simulation
% ********************
caseFolder = [sgvc '_res=' num2str(rez) '_size=' num2str(xArea) 'x' num2str(yArea) 'x' num2str(zArea) '_Lb=' num2str(period) '_d=' num2str(depth) '_w=' num2str(width)];
runMeep = input('run: start Meep simulation? y/n [y] ','s');
if (isempty(runMeep) || runMeep(1) == 'Y' || runMeep(1) == 'y')
    meepControlScript
    system(['mkdir -v -p ' caseFolder]);
    system(['mv ' sgvc '.ctl ' caseFolder]);
    system(['export GUILE_WARN_DEPRECATED="no" && cd ' caseFolder ' && time mpirun -np 4 /usr/bin/meep-mpi-default ' sgvc '.ctl']);
end

%% Import Data, create meshgrid & matrices
% ****************************************
close 'all'
epsZYX = h5read([caseFolder '/' sgvc '-eps-000000.00.h5'],'/eps');
epsYXZ = permute(epsZYX,[2 3 1]);
xMesh = size(epsYXZ,2);     % =xnpts
yMesh = size(epsYXZ,1);     % =ynpts
zMesh = size(epsYXZ,3);     % =znpts
[meshGridX,meshGridY] = meshgrid(1:xMesh,1:yMesh);
meshGridTH = atan2(meshGridY,meshGridX);
meshGridR = sqrt(meshGridX.^2+meshGridY.^2);
for i = 1:size(zOutput,2)
    if force_complex
        Ex(:,:,i) = permute(h5read([caseFolder '/' sgvc '-Ex_z=' num2str(zOutput(i)) '.h5'],'/ex.r') ...
            + h5read([caseFolder '/' sgvc '-Ex_z=' num2str(zOutput(i)) '.h5'],'/ex.i').*1j,[2 3 1]);
        Ey(:,:,i) = permute(h5read([caseFolder '/' sgvc '-Ey_z=' num2str(zOutput(i)) '.h5'],'/ey.r') ...
            + h5read([caseFolder '/' sgvc '-Ey_z=' num2str(zOutput(i)) '.h5'],'/ey.i').*1j,[2 3 1]);
    else
        Ex(:,:,i) = permute(h5read([caseFolder '/' sgvc '-Ex_z=' num2str(zOutput(i)) '.h5'],'/ex'),[3 2 1]);
        Ey(:,:,i) = permute(h5read([caseFolder '/' sgvc '-Ey_z=' num2str(zOutput(i)) '.h5'],'/ey'),[3 2 1]);        
    end
    Ete(:,:,i) = Ex(:,:,i).*-sin(lp/2*meshGridTH) + Ey(:,:,i).*cos(lp/2*meshGridTH);
    Etm(:,:,i) = Ex(:,:,i).* cos(lp/2*meshGridTH) + Ey(:,:,i).*sin(lp/2*meshGridTH);
end
phix = mod(angle(Ex),2*pi);
phiy = mod(angle(Ey),2*pi);
phitetm = mod(angle(Ete)-angle(Etm)-pi/2,2*pi);
intens = abs(Ex).^2+abs(Ey).^2;



%% Figures: structure
% *******************
close 'all'
xMid = xArea/2*rez+1;
yMid = yArea/2*rez+1;
zMid = (zArea/2+zCenter)*rez;
zStart = (zArea/2+zCenter+depth/2)*rez-0.5;
zEnd = (zArea/2+zCenter-depth/2)*rez-0.5;
newSurf
title('structure x-slice')
hS = surfc((epsZYX(:,:,round(yMid))-1)./(nmat^2-1));
%set(hS,'edgeColor','none')
%shading interp
caxis([0 1])
axis([1 yMesh 1 zMesh])
print('-depsc2',[caseFolder '/' sgvc '_epsYZ.eps'],'-r300');
newSurf
title('structure z-slice')
hS = surfc((epsYXZ(:,:,round(zStart))-1)./(nmat^2-1));
caxis([0 1])
axis([1 xMesh 1 yMesh])
print('-depsc2',[caseFolder '/' sgvc '_epsXY.eps'],'-r300');
newSurf
title('\Phi_{TE-TM}')
hS = surfc(phitetm(:,:,1)./pi);
caxis([0 2])
axis([1 xMesh 1 yMesh])
print('-depsc2',[caseFolder '/' sgvc '_phiTETM.eps'],'-r300');
newSurf
title('\Phi_{TE-TM} Residuals')
hS = surfc(abs(phitetm(:,:,1)-pi)./pi);
caxis([0 .5])
axis([1 xMesh 1 yMesh])
print('-depsc2',[caseFolder '/' sgvc '_phiTETMrez.eps'],'-r300');
newSurf
title('\Phi_{TE-TM} Residuals (scaled)')
hS = surfc(abs(phitetm(:,:,1)-pi)./pi);
caxis([0 .2])
axis([1 xMesh 1 yMesh])
print('-depsc2',[caseFolder '/' sgvc '_phiTETMrezScaled.eps'],'-r300');

% drawPhiPan = input('fig: draw phase ramp? y/n [y] ','s');
% if isempty(drawPhiPan) || drawPhiPan(1) == 'Y' || drawPhiPan(1) == 'y'
%     
% end

% %% Figures
% % ********
% makePNG = input('fig: create PNG figures? y/n [y] ','s');
% if isempty(makePNG) || makePNG(1) == 'Y' || makePNG(1) == 'y'
%     system(['cd ' caseFolder ' && h5topng ' sgvc '-eps*.h5 -o ' sgvc '-eps-xMid_x=' num2str(xMid) '.png -x ' num2str(xMid) ' -c bluered']);
%     system(['cd ' caseFolder ' && h5topng ' sgvc '-eps*.h5 -o ' sgvc '-eps-zEnd_z=' num2str(zEnd) '.png -z ' num2str(zEnd) ' -c bluered']);
%     system(['cd ' caseFolder ' && h5topng ' sgvc '-eps*.h5 -o ' sgvc '-eps-zMid_z=' num2str(zMid) '.png -z ' num2str(zMid) ' -c bluered']);
%     system(['cd ' caseFolder ' && h5topng ' sgvc '-eps*.h5 -o ' sgvc '-eps-zStart_z=' num2str(zStart) '.png -z ' num2str(zStart) ' -c bluered']);
% 	system(['cd ' caseFolder ' && h5topng ' sgvc '-E*.h5 -z 0 -c jet']);
% end

