%% FDTD examples with MEEP
% by Lorenzo 2020 based on Christian SGVC simulator 2014
% this script performs some example FDTD simulations with MEEP
% ( abbr 'eps' may be used as dielectric permittivity or eps file ext )

% currently only case 'hwp', 'fqpm' updated -> extend to all cases (new
%   architecture?)


%% setup and case select

clear all, close all
maxIter=100;

% choose example - see description in following section
expl='hwp';

% expls 'waveguide' and 'bent-waveguide' from meep tutorial at
%   https://meep.readthedocs.io/en/latest/Scheme_Tutorials/Basics/
% expl 'pw' from pw-source.ctl at meep-github


%% create meep ctl file

% automated creation of ctl file not yet implemented, create ctl file by 
%   hand and place it in 'caseFolder'; EXCEPT case 'hwp'
% name of 'caseFolder' is arbitrary

rez=10; % in all cases except hwp (changed later)
switch expl
    case 'nth'
        caseFolder = 'NTH_res=10_size=16x8xno_Lb=2orso';
        % no structure - play with ctl file adding and removing ptlike or
        %   extended sources
        % atm 2d, source polarized in z direction to create plane wave
        % records epsilon at beginning, Ez at end
    case 'waveguide'
        caseFolder = 'WG_res=10_size=16x8xno_Lb=2orso_d=';
        % rectangular waveguide in 2d (Glasfaser) of size (infx1xinf) with
        %   pointlike source at left side 1um from border
        % records epsilon at beginning, Ez at end
    case 'bent-waveguide'
        caseFolder = 'BWG_res=10_size=16x16xno_Lb=2orso_d=';
        % rectangular waveguide in 2d with 90d bent at (3.5,-3.5) with
        %   extended source at left side 1um from border turned on slowly
        % records epsilon at beginning, Ez every .6 (timeunits)
    case 'pw'
        caseFolder = 'PW_res=10';
        % plane wave in 2d coming from bottom left corner made of 2 
        %   extended sources at whole left and whole bottom side - meep
        %   example
        % records epsilon at beginning, Ez at end
    case 'grating' % change rez in ctl scritpt and wavelength
        caseFolder = 'G_res=30';
        % no structure (i.e. no grating) - extended rectangular source (2d)
        %   in 3d at top z layer resulting in plane wave propagating 
        %   downwards polarized in y direction - with complex fields, i.e.
        %   Ey = (re,im) ; no x component ; takes ca 20min for rez=30
        % records epsilon at beginning, Ey at end
    case 'grating-on'
        caseFolder = 'GO_res=30';
        % grating in y direction at z=o plane with source of 'grating' case
        %   + x component at -lam/4 downwards -> LHC wave
        % records epsilon at beginning, Ey at end
    case 'hwp'
        % half wave plate, i.e. SG for L band with params Lb=lam/n=1.47,
        %   lw=.65, d=5.50, sw=0
        % source lam=3.5 no bandwidth (center of L band) pw 45deg pol (no
        %   LHC or RHC)
        % --- set parameters
        rez=16;
        lam=3.5;
        xsize=12;
        ysize=12;
        zsize=14;
        Lb=1.21; % must be \lt lam/n, n=2.38xx to meet SG condition
        nLb=ceil(xsize./Lb); % round up
        lw=.7876;
        d=4.5753;
        runtime=72;
        isintegrated='true'; % necessary for pw source extending into pml
        comment='source=45both_isintegrated=true';
        caseFolder = ['HWP_res=' num2str(rez) '_lam=' num2str(lam) '_size=' num2str(xsize) 'x' num2str(ysize) 'x' num2str(zsize) '_Lb=' num2str(Lb) '_lw=' num2str(lw) '_d=' num2str(d) '_t=' num2str(runtime) '_' comment];
        % --- write ctl file
        f=fopen(['expl_' expl '.ctl'],'w');
        fprintf(f,['(set! resolution ' num2str(rez) ')\n']);
        fprintf(f,'(set! force-complex-fields? true)\n');
        fprintf(f,'(set! pml-layers (list\n   (make pml (thickness 1))\n))\n');
        fprintf(f,'(set! sources (list\n');
        fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')(is-integrated? ' isintegrated ')))(component Ex)(center 0 0 ' num2str(zsize./2-1) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']);
        fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')(is-integrated? ' isintegrated ')))(component Ey)(center 0 0 ' num2str(zsize./2-1-0.*lam./4) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']);
        fprintf(f,'))');
        fprintf(f,['(set! geometry-lattice (make lattice (size ' num2str(xsize) ' ' num2str(ysize) ' ' num2str(zsize) ')))\n']);
        fprintf(f,'(set! geometry (list\n');
        for i=1:nLb
            fprintf(f,['   (make block (center ' num2str(-xsize./2+(i-.5).*Lb) ' 0 ' num2str(zsize./2-3.25-d./2) ')(size ' num2str(lw) ' ' num2str(ysize) ' ' num2str(d) ')(material (make dielectric (index 2.3812))))\n']);
        end % xcenter from period; zcenter is 3.25 space on top, and 'd' from top
        fprintf(f,['   (make block (center 0 0 ' num2str(-(3.25+d)./2) ')(size ' num2str(xsize) ' ' num2str(ysize) ' ' num2str(zsize-(3.25+d)) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,'))\n'); % zcenter is 3.25+d+(1./2 from remaining space)
        fprintf(f,['(run-until ' num2str(runtime) '\n   (at-beginning\n      output-epsilon\n   )\n   (at-end\n      output-efield-x\n      output-efield-y\n   )\n)']);
        fclose(f);
        % --- copy to casefolder
        system(['mkdir -v -p ' caseFolder]);
        system(['mv expl_' expl '.ctl ' caseFolder]);
    case 'fqpm'
        % four quadrant phase mask, i.e. discretized vortez phase ramp for 
        %   lam=3.5, four steps of pi/2 each
        % source lam=3.5 no bandwidth (center of L band) pw 45deg pol (no
        %   LHC or RHC)
        % --- set parameters
        rez=6;
        lam=3.5;
        d=7.60; %2*[1,3,5,...]*lam/2/(n-1)
        xsize=12;
        ysize=12;
        zsize=16;
        zoff=3.25;%6.25;16.25;
        runtime=36;%=(zoff-1)*n_air+(zsize-zoff-1pml)*n_diamond+20%puffer
        isintegrated='false'; % necessary for pw source extending into pml
        comment=['d=' num2str(d) '_zoff=' num2str(zoff) '_source=45both'];
        caseFolder = ['FQPM_res=' num2str(rez) '_lam=' num2str(lam) '_size=' num2str(xsize) 'x' num2str(ysize) 'x' num2str(zsize) '_t=' num2str(runtime) '_' comment '_is-integrated=' isintegrated];
        % write ctl file
        f=fopen(['expl_' expl '.ctl'],'w');
        fprintf(f,['(set! resolution ' num2str(rez) ')\n']);
        fprintf(f,'(set! force-complex-fields? true)\n');
        fprintf(f,'(set! pml-layers (list\n   (make pml (thickness 1))\n))\n');
        fprintf(f,'(set! sources (list\n');
        %fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')))(component Ex)(center 0 0 ' num2str(zsize./2-1) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']);
        %fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')))(component Ey)(center 0 0 ' num2str(zsize./2-1-0.*lam./4) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']);
        fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')(is-integrated? ' isintegrated ')))(component Ex)(center 0 0 ' num2str(zsize./2-1) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']);
        fprintf(f,['   (make source (src (make continuous-src (frequency ' num2str(1./lam) ')(is-integrated? ' isintegrated ')))(component Ey)(center 0 0 ' num2str(zsize./2-1-0.*lam./4) ')(size ' num2str(xsize) ' ' num2str(ysize) ' 0)(amplitude 1))\n']); % for both Ex and Ey!
        fprintf(f,'))');
        fprintf(f,['(set! geometry-lattice (make lattice (size ' num2str(xsize) ' ' num2str(ysize) ' ' num2str(zsize) ')))\n']);
        fprintf(f,'(set! geometry (list\n');
        fprintf(f,['   (make block (center ' num2str(-xsize./4) ' ' num2str(-ysize./4) ' ' num2str(zsize./2-zoff-d+.25.*d./2) ')(size ' num2str(xsize./2) ' ' num2str(ysize./2) ' ' num2str(.25.*d) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,['   (make block (center ' num2str(-xsize./4) ' ' num2str(+ysize./4) ' ' num2str(zsize./2-zoff-d+.50.*d./2) ')(size ' num2str(xsize./2) ' ' num2str(ysize./2) ' ' num2str(.50.*d) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,['   (make block (center ' num2str(+xsize./4) ' ' num2str(+ysize./4) ' ' num2str(zsize./2-zoff-d+.75.*d./2) ')(size ' num2str(xsize./2) ' ' num2str(ysize./2) ' ' num2str(.75.*d) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,['   (make block (center ' num2str(+xsize./4) ' ' num2str(-ysize./4) ' ' num2str(zsize./2-zoff-d+1.0.*d./2) ')(size ' num2str(xsize./2) ' ' num2str(ysize./2) ' ' num2str(1.0.*d) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,['   (make block (center 0 0 ' num2str((-zoff-d)./2) ')(size ' num2str(xsize) ' ' num2str(ysize) ' ' num2str(zsize-zoff-d) ')(material (make dielectric (index 2.3812))))\n']);
        fprintf(f,'))\n');
        fprintf(f,['(run-until ' num2str(runtime) '\n   (at-beginning\n      output-epsilon\n   )\n   (at-end\n      output-efield-x\n      output-efield-y\n   )\n)']);
        fclose(f);
        % copy to casefolder
        system(['mkdir -v -p ' caseFolder]);
        system(['mv expl_' expl '.ctl ' caseFolder]);
end


%% meep simulation

disp('starting meep simulation')
system(['export GUILE_WARN_DEPRECATED="no" && cd ' caseFolder ' && time /usr/local/bin/meep expl_' expl '.ctl']);
disp('finished meep simulation')


%% read resulting h5data

close all

% read epsilon distribution at beginning of simulation 000000.00
epsdata1=h5read([caseFolder '/expl_' expl '-eps-000000.00.h5'],'/eps');
epsdata=epsdata1; %for 2d case. for 3d case epsdata is redefined later

% read E field at different times (depending on case)
% define size xyzMesh of simulation grid (for plotting)
% compute case specific quantities, i.e. amplitude and phase for complex
%   fields
switch expl
    case 'nth'
        ezdata=h5read([caseFolder '/expl_' expl '-ez-000200.00.h5'],'/ez');
        xMesh=160;
        yMesh=160;
    case 'waveguide'
        ezdata=h5read([caseFolder '/expl_' expl '-ez-000200.00.h5'],'/ez');
        xMesh=160;
        yMesh=80;
   case 'bent-waveguide'
        ezdata1=h5read([caseFolder '/expl_' expl '-ez.h5'],'/ez');
        ezdata=permute(ezdata1,[2 3 1]); % 'movie' - '1' is time dimension
        xMesh=160;
        yMesh=160;
   case 'pw'
        ezdata=h5read([caseFolder '/expl_' expl '-ez-000400.00.h5'],'/ez');
        xMesh=130;
        yMesh=130;
   case {'grating','grating-on'} % with complex fields; get amplitude and phase
        epsdata=squeeze(epsdata1(30,:,:)); % redef; (Z,Y,X) get a zsclice
        eyre1=h5read([caseFolder '/expl_' expl '-ey-000024.00.h5'],'/ey.r');
        eyim1=h5read([caseFolder '/expl_' expl '-ey-000024.00.h5'],'/ey.i');
        eydata1=eyre1+1j.*eyim1;
        eyam1=abs(eydata1);
        eyph1=angle(eydata1);
        eyam=squeeze(eyam1(30,:,:));
        eyph=squeeze(eyph1(30,:,:)); % x part only for grating-on
        exre1=h5read([caseFolder '/expl_' expl '-ex-000024.00.h5'],'/ex.r');
        exim1=h5read([caseFolder '/expl_' expl '-ex-000024.00.h5'],'/ex.i');
        exdata1=exre1+1j.*exim1;
        exam1=abs(exdata1);
        exph1=angle(exdata1);
        exam=squeeze(exam1(30,:,:));
        exph=squeeze(exph1(30,:,:));
        xMesh=100;
        yMesh=100; % eig. 80
        rez=30;       % 'grating' with higher rez
        xMesh=3.*100; %
        yMesh=3.*80;  %
   case {'hwp','fqpm'} % with complex fields; get amplitude and phase
        epsdata=squeeze(epsdata1(30,:,:)); % redef; (Z,Y,X) get a zsclice
        eyre1=h5read([caseFolder '/expl_' expl '-ey-' num2str(runtime,'%06.f') '.00.h5'],'/ey.r');
        eyim1=h5read([caseFolder '/expl_' expl '-ey-' num2str(runtime,'%06.f') '.00.h5'],'/ey.i');
        eydata1=eyre1+1j.*eyim1;
        eyam1=abs(eydata1);
        eyph1=angle(eydata1);
        eyam=squeeze(eyam1(30,:,:));
        eyph=squeeze(eyph1(30,:,:)); % x part only for grating-on
        exre1=h5read([caseFolder '/expl_' expl '-ex-' num2str(runtime,'%06.f') '.00.h5'],'/ex.r');
        exim1=h5read([caseFolder '/expl_' expl '-ex-' num2str(runtime,'%06.f') '.00.h5'],'/ex.i');
        exdata1=exre1+1j.*exim1;
        exam1=abs(exdata1);
        exph1=angle(exdata1);
        exam=squeeze(exam1(30,:,:));
        exph=squeeze(exph1(30,:,:));
        dphi1=mod(eyph1-exph1,2.*pi);
        xMesh=xsize.*rez;
        yMesh=ysize.*rez;
        q=(eyam1./exam1).^2;
        e=dphi1-pi;
        N=((1-sqrt(q)).^2+e.^2.*sqrt(q))./(1+sqrt(q)).^2; %Mawet05
end


%% plot structure

% plot epsilon distribution - plot a cut at z=sth in 3d case
switch expl
    case {'hwp','fqpm'}
        % already plotted within gif later
    otherwise
        vis='on';
        mnewsurf % based on newsurf by Christian
        title('epsilon distribution (structure)')
        colormap(flipud(gray(256)))
        caxis([0 12]) % epsilon = 12 (medium) and 1 (air) are used
        grid off
        surf(epsdata,'edgecolor','none','FaceAlpha',0.33) 
        % to avoid drawing pixel borders add (...,'edgecolor','none')
        print('-depsc2',[caseFolder '/expl_' expl '_eps.eps'],'-r300');
end

% plot E field (case specific)
disp('plotting E field...');
switch expl
    % === Half Wave Plate - save several z slices as gif
    % === FQPM - comments do not hold
    case {'hwp','fqpm'}
        sca=1.5; % convenient color bar scale for amplitude plots
        vis='off';
        isfirstloop=1;
        for i=1:1:rez*zsize
            disp([num2str(i) ' von ' num2str(rez*zsize)]);
            % --- plot x intensity - should be 1 (const) everywhere
            mnewsurf
            title(['Ex intensity - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 sca])
            colormap(jet)
            surf(squeeze(exam1(rez.*zsize-i+1,:,:)));
            frame1=getframe(h);
            % --- do not plot y intensity
            % --- plot q=TE/TM=I_y/I_x; y is parallel to grooves
            mnewsurf
            title(['q=I_y/I_x=TE/TM - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 1.5])
            colormap(jet)
            surf(squeeze(q(rez.*zsize-i+1,:,:)));
            frame2=getframe(h);
            % --- plot epsilon - should be grating structure
            mnewsurf
            title(['epsilon dist - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([1 6])
            colormap(flipud(gray))
            surf(squeeze(epsdata1(rez.*zsize-i+1,:,:)));
            frame3=getframe(h);
            % --- plot x phase - should go from 0 to 2pi continuously
            mnewsurf
            %phixmean=mean(mod(squeeze(exph1(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),2.*pi),'all');
            phixmedian=median(mod(squeeze(exph1(rez*zsize-i+1,:,:)),2.*pi),'all');
            title(['Ex phase (median) = ' num2str(phixmedian) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 2.*pi])    %^median uses whole xyrange, not cut by 15 (as was the case with mean)
            colormap default
            colorbar('Ticks',[0,pi./2,pi,3.*pi./2,2.*pi],'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            surf(mod(squeeze(exph1(rez.*zsize-i+1,:,:)),2.*pi));
            frame4=getframe(h);
            % --- do not plot y phase
            % --- plot e=phase diff minus pi - should be zero after SG
            mnewsurf
            %dphimean=mean(squeeze(dphi1(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),'all');
            dphimedian=median(squeeze(dphi1(rez*zsize-i+1,:,:)),'all');
            title(['e = dphi-pi (median) = ' num2str(dphimedian-pi) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([-1.0.*pi 1.0.*pi])
            colormap(redblue)
            colorbar('Ticks',[-pi,-pi./2,0,pi./2,pi],'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
            surf(squeeze(e(rez.*zsize-i+1,:,:)));
            frame5=getframe(h);
            % --- plot Null depth N=complicated formula with q and e
            mnewsurf
            %Nmean=mean(squeeze(N(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),'all');
            Nmedian=median(squeeze(N(rez*zsize-i+1,:,:)),'all');
            title(['Null depth (median) = ' num2str(Nmedian) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 1])
            colormap default
            surf(squeeze(N(rez.*zsize-i+1,:,:)));
            frame6=getframe(h);
            % --- merge 6 frames to image
            im1=frame2im(frame1);
            im2=frame2im(frame2);
            im3=frame2im(frame3);
            im4=frame2im(frame4);
            im5=frame2im(frame5);
            im6=frame2im(frame6);
            imm=montage({im1,im2,im3,im4,im5,im6});
            im=imm.CData;
            % --- append to gif
            [imind,cm] = rgb2ind(im,256); 
            if isfirstloop == 1
                imwrite(imind,cm,[caseFolder '/data_as_gif_for_diff_zslice.gif'],'gif', 'Loopcount',inf);
                isfirstloop=0;
            else 
                imwrite(imind,cm,[caseFolder '/data_as_gif_for_diff_zslice.gif'],'gif','WriteMode','append');
            end
        close all
        end
    case {'hwp2','fqpm2'} % backup of previous case with 6 images in gif
        sca=1.5; % convenient color bar scale for amplitude plots
        vis='off';
        isfirstloop=1;
        for i=1:1:rez*zsize
            disp([num2str(i) ' von ' num2str(rez*zsize)]);
            % --- plot intensity - should be 1 (const) everywhere
            mnewsurf
            title(['Ex intensity - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 sca])
            colormap(jet)
            surf(squeeze(exam1(rez.*zsize-i+1,:,:)));
            frame1=getframe(h);
            % --- plot intensity - should be 1 (const) everywhere
            mnewsurf
            title(['Ey intensity - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 sca])
            colormap(jet)
            surf(squeeze(eyam1(rez.*zsize-i+1,:,:)));
            frame2=getframe(h);
            % --- plot epsilon - should be grating structure
            mnewsurf
            title(['epsilon dist - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([1 6])
            colormap(flipud(gray))
            surf(squeeze(epsdata1(rez.*zsize-i+1,:,:)));
            frame3=getframe(h);
            % --- plot phase - should go from 0 to 2pi continuously
            mnewsurf
            phixmean=mean(mod(squeeze(exph1(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),2.*pi),'all');
            title(['Ex phase = ' num2str(phixmean) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 2.*pi])
            colormap default
            colorbar('Ticks',[0,pi./2,pi,3.*pi./2,2.*pi],'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            surf(mod(squeeze(exph1(rez.*zsize-i+1,:,:)),2.*pi));
            frame4=getframe(h);
            % --- plot phase - should go from 0 to 2pi continuously
            mnewsurf
            phiymean=mean(mod(squeeze(eyph1(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),2.*pi),'all');
            title(['Ey phase = ' num2str(phiymean) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 2.*pi])
            colormap default
            colorbar('Ticks',[0,pi./2,pi,3.*pi./2,2.*pi],'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            surf(mod(squeeze(eyph1(rez.*zsize-i+1,:,:)),2.*pi));
            frame5=getframe(h);
            % --- plot phase diff - should go from 0 to 2pi continuously
            mnewsurf
            dphimean=mean(squeeze(dphi1(rez*zsize-i+1,1.5*rez:(ysize-1.5)*rez,1.5*rez:(xsize-1.5)*rez)),'all');
            title(['Delta phase = ' num2str(dphimean) ' - z=' num2str((rez.*zsize-i+1)./rez-zsize./2)])
            caxis([0 2.*pi])
            colormap default
            colorbar('Ticks',[0,pi./2,pi,3.*pi./2,2.*pi],'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            surf(squeeze(dphi1(rez.*zsize-i+1,:,:)));
            frame6=getframe(h);
            % --- merge 6 frames to image
            im1=frame2im(frame1);
            im2=frame2im(frame2);
            im3=frame2im(frame3);
            im4=frame2im(frame4);
            im5=frame2im(frame5);
            im6=frame2im(frame6);
            imm=montage({im1,im2,im3,im4,im5,im6});
            im=imm.CData;
            % --- append to gif
            [imind,cm] = rgb2ind(im,256); 
            if isfirstloop == 1
                imwrite(imind,cm,[caseFolder '/data_as_gif_for_diff_zslice.gif'],'gif', 'Loopcount',inf);
                isfirstloop=0;
            else 
                imwrite(imind,cm,[caseFolder '/data_as_gif_for_diff_zslice.gif'],'gif','WriteMode','append');
            end
        close all
        end
    case {'grating','grating-on'} % get several z slices as gif
        % some 5.*() for 'grating' with higher rez
        sca=max(abs(eyam1(:)))./2; % arbitrary 
        vis='off';
        isfirstloop=1;
        for i=5:2:110
            mnewsurf % plot intensity - should be 1 everywhere
            title(['Ey intensity - z=' num2str(3.*(120-i))])
            caxis([0 sca])
            colormap(jet)
            surf(squeeze(eyam1(3.*(120-i),:,:)));
            frame1=getframe(h);
            mnewsurf % plot intensity - should be 1 everywhere
            title(['epsilon dist - z=' num2str(3.*(120-i))])
            caxis([1 6])
            colormap(flipud(gray))
            surf(squeeze(epsdata1(3.*(120-i),:,:)));
            frame2=getframe(h);
            mnewsurf % plot phase - should go from 0 to 2pi continuously
            title(['Ey phase     - z=' num2str(3.*(120-i))])
            caxis([0 2.*pi])
            colormap default
            surf(mod(squeeze(eyph1(3.*(120-i),:,:)),2.*pi)); % to avoid drawing pixel borders add (...,'edgecolor','none')
            frame3=getframe(h);
            %mnewsurf % plot 3d structure with current z slice cut
            %title('epsilon ')             % 4th image could be 3d struc-
            %caxis([0 2.*pi])              %   ture, but quite complicated
            %colormap default              %   and costly in comput. time
            %surf(mod(squeeze(eyph1(3.*(120-i),:,:)),2.*pi)); % to avoid drawing pixel borders add (...,'edgecolor','none')
            %frame4=getframe(h);           %
            im1=frame2im(frame1);
            im2=frame2im(frame2);
            im3=frame2im(frame3);
            imm=montage({im1,im2,im3});
            im=imm.CData;
            [imind,cm] = rgb2ind(im,256); 
            if isfirstloop == 1
                imwrite(imind,cm,[caseFolder '/yamplitudeasgifdiffzslice.gif'],'gif', 'Loopcount',inf); 
            else 
                imwrite(imind,cm,[caseFolder '/yamplitudeasgifdiffzslice.gif'],'gif','WriteMode','append'); 
            end
            %mnewsurf % plot phase - should go from 0 to 2pi continuously
            %title(['Ey phase - z=' num2str(3.*(120-i))])
            %caxis([0 2.*pi])
            %colormap default
            %surf(mod(squeeze(eyph1(3.*(120-i),:,:)),2.*pi)); % to avoid drawing pixel borders add (...,'edgecolor','none')
            %frame=getframe(h);
            %im=frame2im(frame);
            %[imind,cm] = rgb2ind(im,256); 
            %if isfirstloop == 1
            %    imwrite(imind,cm,[caseFolder '/yphaseasgifdiffzslice.gif'],'gif', 'Loopcount',inf); 
            %else 
            %    imwrite(imind,cm,[caseFolder '/yphaseasgifdiffzslice.gif'],'gif','WriteMode','append'); 
            %end
            isfirstloop=0;
        close all
        end
    case 'nth'
        sca=max(abs(ezdata(:)))./3.;
        mnewsurf
        title('Ez field')
        caxis([-sca sca])
        colormap(redblue)
        surf(ezdata) % to avoid drawing pixel borders add (...,'edgecolor','none')
        print('-depsc2',[caseFolder '/expl_' expl '_ez.eps'],'-r300');
    case 'pw'
        sca=max(abs(ezdata(:)))./3.;
        mnewsurf
        title('Ez field')
        caxis([-sca sca])
        colormap(redblue)
        surf(ezdata) % to avoid drawing pixel borders add (...,'edgecolor','none')
        print('-depsc2',[caseFolder '/expl_' expl '_ez.eps'],'-r300');
    case 'waveguide' % 2d
        mnewsurf
        title('Ez field')
        colormap(redblue)
        surf(ezdata) % to avoid drawing pixel borders add (...,'edgecolor','none')
        print('-depsc2',[caseFolder '/expl_' expl '_ez.eps'],'-r300');
    case 'bent-waveguide' % 2d get time series as gif
        sca=max(abs(ezdata(:)))./10;
        vis='off';
        isfirstloop=1;
        for i=150:1:166
            mnewsurf
            title('Ez field')
            caxis([-sca sca])
            colormap(redblue)
            surf(ezdata(:,:,i)); % to avoid drawing pixel borders add (...,'edgecolor','none')
            print('-depsc2',[caseFolder '/expl_' expl '_ez.eps'],'-r300');
            frame=getframe(h);
            im=frame2im(frame);
            [imind,cm] = rgb2ind(im,256); 
            if isfirstloop == 1
                imwrite(imind,cm,[caseFolder '/zdataasgif.gif'],'gif', 'Loopcount',inf); 
            else 
                imwrite(imind,cm,[caseFolder '/zdataasgif.gif'],'gif','WriteMode','append'); 
            end
            isfirstloop=0;
        close all
        end
end

disp('... E field done.')