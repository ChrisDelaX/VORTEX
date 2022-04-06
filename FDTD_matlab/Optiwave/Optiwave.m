clear all
tic

% Input parameters - import data
% ------------------------------
sgvc = 'HWP_horiz_.05-.05-500-5_Lb1.4_d5.15_w.9/output';%'SGVC4_disc_2048x30';
LP = 0;%4;%2;%
thHWP = 1*pi/2;  %vertical=0 ; horizontal = pi/2
meshSize = 0.05;%0.1;%0.2;%
xMeshData = 601;%;%391;%451;%241;%181;%361;%721;%257;%513;%1025;%451;%
yMeshData = xMeshData;
xMirror = 1;%0;% 
yMirror = 1;%0;% 
xyMesh = 601;%601;%391;%451;%1025;%513;1025;%
figMesh = 151;%    %(figMesh-1) MUST be a fraction of (xyMesh-1), and < 175 !!!
ringMesh = 15;%32;%13;%30;%32;%25;%   % fraction of (xyMesh-1)/2
%quadrants = 4;% 0;%
Extmp = fscanf(fopen(sprintf('%s/Ex.txt',sgvc)), '%g %g', [2 inf])';
Eytmp = fscanf(fopen(sprintf('%s/Ey.txt',sgvc)), '%g %g', [2 inf])';
Eztmp = fscanf(fopen(sprintf('%s/Ez.txt',sgvc)), '%g %g', [2 inf])';
if size(Extmp,1) ~= xMeshData*yMeshData
    error('Wrong number of mesh units.')
end
xMesh = (xyMesh-1)/(1+yMirror)+1;
yMesh = (xyMesh-1)/(1+xMirror)+1;
prealloc % Preallocate matrixes
for ii=1:yMesh
    for jj=1:xMesh
        kk = yMesh-ii+1;  % Invert y-axis (colums) because of Matlab
        ll = jj;%xMesh-jj+1;  % Invert x-axis (colums) because of Optiwave -- ok 
        linePos = (ii+yMeshData-yMesh-1)*xMeshData+jj+yMeshData-yMesh;%+xMeshData-xMesh;
        Ex(kk,ll) = Extmp(linePos,1)+1i*Extmp(linePos,2);
        Ey(kk,ll) = Eytmp(linePos,1)+1i*Eytmp(linePos,2);
        Ez(kk,ll) = Eztmp(linePos,1)+1i*Eztmp(linePos,2);
    end
end

%    Ex = flipdim(Ex,2);  %rot90 %flipdim
%    Ey = flipdim(Ey,2);  %rot90 %flipdim
%    Ez = flipdim(Ez,2);  %rot90 %flipdim
if yMirror == 1 || xMirror == 1
    if yMirror == 0  % only xMirror == 1
        Ex(end:2*end-1,:) = rot90((Ex),2);
        Ey(end:2*end-1,:) = rot90((Ey),2);
        Ez(end:2*end-1,:) = rot90((Ez),2);
    elseif xMirror == 0  % only yMirror == 1
        Ex(:,end:2*end-1) = rot90((Ex),2);
        Ey(:,end:2*end-1) = rot90((Ey),2);
        Ez(:,end:2*end-1) = rot90((Ez),2);
    else % both == 0
        Ex(end:2*end-1,:) = rot90((Ex),1);
        Ey(end:2*end-1,:) = rot90((Ey),1);
        Ez(end:2*end-1,:) = rot90((Ez),1);
        Ex(:,end:2*end-1) = rot90((Ex),2);
        Ey(:,end:2*end-1) = rot90((Ey),2);
        Ez(:,end:2*end-1) = rot90((Ez),2);
    end
end
        
        %

%         
%     else if xMirror == 0
%         else
% else if end
    
% % for AFTA
% meshSize = 0.02;%
% xyMesh = 641;%
% figMesh = 161;%  %(figMesh-1) MUST be a fraction of (xyMesh-1), and < 175 !!!
% ringMesh = 32;%25;%   % fraction of (xyMesh-1)/2


% Create meshGrid and matrices
% ----------------------------
center = (xyMesh-1)/2+1;
ringWidth = round((center-1)/ringMesh);
[meshGridX,meshGridY] = meshgrid(-center+1:center-1);
[meshGridTH,meshGridR] = cart2pol(meshGridX,meshGridY);
realGrid = linspace(-meshSize*(center-1),meshSize*(center-1),xyMesh);
realGridFig = linspace(-meshSize*(center-1),meshSize*(center-1),figMesh);
for ii=1:xyMesh
    for jj=1:xyMesh
        Ete(ii,jj) = Ex(ii,jj)*-sin(LP/2*meshGridTH(ii,jj)+thHWP)+Ey(ii,jj)*cos(LP/2*meshGridTH(ii,jj)+thHWP);
        Etm(ii,jj) = Ex(ii,jj)* cos(LP/2*meshGridTH(ii,jj)+thHWP)+Ey(ii,jj)*sin(LP/2*meshGridTH(ii,jj)+thHWP);
        phix(ii,jj) = mod(angle(Ex(ii,jj))+2*pi,2*pi);
        phiy(ii,jj) = mod(angle(Ey(ii,jj))+2*pi,2*pi);
        phiz(ii,jj) = mod(angle(Ez(ii,jj))+2*pi,2*pi);
        phite(ii,jj) = mod(angle(Ete(ii,jj))+2*pi,2*pi);
        phitm(ii,jj) = mod(angle(Etm(ii,jj))+2*pi,2*pi);
        phitetm(ii,jj) = mod(phite(ii,jj)-phitm(ii,jj)+2*pi,2*pi)-pi/2;
        Ax(ii,jj) = abs(Ex(ii,jj));
        Ay(ii,jj) = abs(Ey(ii,jj));
        Az(ii,jj) = abs(Ez(ii,jj));
        Int(ii,jj) = Ax(ii,jj)^2+Ay(ii,jj)^2+Az(ii,jj)^2;
    end
end

% Reduced matrices for figures
xIncrement = (xyMesh-1)/(figMesh-1);
yIncrement = (xyMesh-1)/(figMesh-1);
iiFig=1;
for ii=1:xIncrement:xyMesh
    jjFig=1;
    for jj=1:yIncrement:xyMesh
        ExFig(iiFig,jjFig) = Ex(ii,jj);
        EyFig(iiFig,jjFig) = Ey(ii,jj);
        EzFig(iiFig,jjFig) = Ez(ii,jj);
        EteFig(iiFig,jjFig) = Ete(ii,jj);
        EtmFig(iiFig,jjFig) = Etm(ii,jj);
        phixFig(iiFig,jjFig) = phix(ii,jj);
        phiyFig(iiFig,jjFig) = phiy(ii,jj);
        phizFig(iiFig,jjFig) = phiz(ii,jj);
        phitetmFig(iiFig,jjFig) = phitetm(ii,jj);
        AxFig(iiFig,jjFig) = Ax(ii,jj);
        AyFig(iiFig,jjFig) = Ay(ii,jj);
        AzFig(iiFig,jjFig) = Az(ii,jj);
        IntFig(iiFig,jjFig) = Int(ii,jj);
        jjFig=jjFig+1;
    end
    iiFig=iiFig+1;
end

% Mean radial (off-axis) intensity
count(1) = 1;
sumInt(1) = Int(center,center);
meanInt(1) = Int(center,center);
offAxis(1)=0;
for ii=1:ringMesh
    count(ii+1) = 0;
    sumInt(ii+1) = 0;
    offAxis(ii+1)=(ii-1/2)*meshSize*ringWidth;
    rmin = (ii-1)*ringWidth;
    rmax = ii*ringWidth;
    for jj=center-ii*ringWidth:center+ii*ringWidth
        for kk=center-ii*ringWidth:center+ii*ringWidth
            if rmin < meshGridR(jj,kk) && meshGridR(jj,kk) <= rmax
                count(ii+1) = count(ii+1)+1;
                sumInt(ii+1) = sumInt(ii+1)+Int(jj,kk);
            end
        end
    end
    meanInt(ii+1) = sumInt(ii+1)/count(ii+1);
end

% Pancharatnam phase
U = 1/sqrt(2)*[1 1i;1 -1i]; % helical-basis transformation matix = [R;L]
%Vrefx = Ex(center,round((xyMesh-1)*1/4+1));  % attention: lign=y , column=x
Vrefx = Ex(center,1);  % attention: lign=y , column=x
Vrefy = Ey(center,1);
Vref = U*[Vrefx;Vrefy];
for ii=1:xyMesh
    for jj=1:xyMesh
        Vpt = U*[Ex(ii,jj);Ey(ii,jj)];
        prodscal = Vref'*Vpt;
        phiPan(ii,jj) = mod(angle(prodscal)+2*pi,2*pi);
        orient(ii,jj) = mod(meshGridTH(ii,jj).*(LP/2),2*pi);
    end
end
iiFig=1;
for ii=1:xIncrement:xyMesh
    jjFig=1;
    for jj=1:yIncrement:xyMesh
        phiPanFig(iiFig,jjFig) = phiPan(ii,jj);
        orientFig(iiFig,jjFig) = orient(ii,jj);
        jjFig=jjFig+1;
    end
    iiFig=iiFig+1;
end

toc
%save(sprintf('%s/data.mat',sgvc))
%save(sprintf('%s/data.fits',sgvc))



% Ax,y,z , Phix,y,z , Phi_TE-TM , PhiPan
% --------------------------------------
FiguresPhases

% Ex,y,z (Re and Im)
% ------------------
FiguresExyz

% Int, log(Int)
% -------------
%FiguresInt

