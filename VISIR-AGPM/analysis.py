# -*- coding: utf-8 -*-
import Tools.ImageProcessing as ImPro
from Tools.ObservingBlock import ObservingBlock
import os.path
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt



""" Input Data """

instruPath = os.path.expandvars('$HOME/INSTRUMENTS/VISIR')
nbin = 1.           # binning parameter for the radial profile
#asecpmm = 1.9484   # cf. Eric Pantin
asecpp = 45.3e-3    # cf. Eric Pantin (arcsec per pixel = asecpmm*.8/34)
LS = 0.9            # Lyot-stop
diam = 8.2          # telescope diameter
chopnod = 178       # chopping-nodding offset (in pixels)
rim = 50            # image region of interest
rbg = 250           # background region of interest

targetNum = 1       # select the set of images corresponding to a run/target
offaxisFit = False  # to find the off-axis PSF centers


""" Targets and filters """

offsets = []
if targetNum == 0:
    target = '2015-03-09_off-axis'
    seq = ['0020','0021','0023','0022','0026','0027','0030','0031','0033','0032','0034',\
            '0035','0037','0036','0038','0039','0041','0040','0042','0043','0044','0045',\
            '0046','0047','0048','0049','0050','0052','0053','0054'] # from Eric
    filter = ['BB-AGPM']*len(seq)
    lam = np.array([12.4]*len(seq))*1e-6
    seeing = np.zeros(np.size(seq))
    # centers: to print offsets as string: ('%.2f' % offsets[i]).zfill(6)
    offsets = np.array([114.85,114.85,114.75,114.75,114.65,114.65,114.55,114.55,114.50,\
            114.50,114.50,114.50,114.50,114.50,114.45,114.45,114.40,114.40,114.35,114.35,\
            114.30,114.30,114.25,114.25,114.20,114.20,114.15,114.15,114.05,114.05]) # from Eric
    pxmm =  (407.218-372.885)/0.8 # pixels per mm offset (407.2-373.9)/0.8 #
    cx = np.array(407.218 - pxmm*(114.85 - offsets), dtype=int) # cf. fit of offsets
    cy = np.array([164]*np.size(seq), dtype=int)
    # specific to this run
    chopnod = 0
    rim = 15
elif targetNum == 1:
    target = 'HD125932'
    seq = ['0009','0010','0011','0012','0013','0014','0015','0016','0017','0018','0019','0020']
    filter = ['BB-AGPM','BB-AGPM','PAH1','PAH1','NB-4QPM1','NB-4QPM1','NeII','NeII',\
            'PAH2','PAH2','NB-4QPM2','NB-4QPM2']
    lam = np.array([12.4,12.4,8.6,8.6,10.65,10.65,12.4,12.4,11.3,11.3,11.3,11.3])*1e-6
    seeing = np.array([1.07,1.07,1.21,1.21,1.16,1.16,1.70,1.70,1.18,1.18,0.84,0.84])
    # centers
    cx = np.array([391]*np.size(seq), dtype=int)
    cy = np.array([437, 437, 549, 549, 437, 437, 549, 549, 549, 549, 437, 437], dtype=int)
elif targetNum == 2:
    target = 'HD139063'
    seq = ['0021','0022','0023','0024','0025','0026','0027','0028','0029','0030']
    filter = ['BB-AGPM','BB-AGPM','PAH1','PAH1','NB-4QPM1','NB-4QPM1','NB-4QPM2',\
            'NB-4QPM2','BB-AGPMbis','BB-AGPMbis']
    lam = np.array([12.4,12.4,8.6,8.6,10.65,10.65,11.3,11.3,12.4,12.4])*1e-6
    seeing = np.array([0.89,0.89,0.90,0.90,0.72,0.72,0.83,0.83,0.96,0.96])
    bws = np.array([1.6,1.6,1.,1.,0,0,0,0,1.6,1.6])*1e-6 # ????
    # centers
    cx = np.array([391]*np.size(seq), dtype=int)
    cy = np.array([437, 437, 549, 549, 437, 437, 437, 437, 437, 437], dtype=int)
elif targetNum == 3:
    target = 'HR6378'
    seq = ['0031','0032']
    filter = ['BB-AGPM','BB-AGPM']
    lam = np.array([12.4,12.4])*1e-6
    seeing = np.array([1.22,1.22])
    # centers
    cx = np.array([391, 391], dtype=int)
    cy = np.array([437, 437], dtype=int)


""" create the Observing Block """

path = os.path.join(instruPath, target)
OB = ObservingBlock(path, seq=seq)
OB.filter = np.array(filter, ndmin=1)
OB.lam = np.array(lam, ndmin=1)
OB.lamD = OB.lam/(diam*LS)*(60*60*180/np.pi)/asecpp     # lam/D in pixel
OB.FWHM = 1*OB.lamD                                     # Full Width Half Max
OB.seeing = np.array(seeing, ndmin=1)
OB.cx, OB.cy = cx, cy
del(cx, cy)
OB.cx2, OB.cy2 = [],[]


""" Analysis """

psfAB = np.empty(((OB.nfiles/2,2*rim+1,2*rim+1)))
psfAB_AGPM = np.empty(((OB.nfiles/2,2*rim+1,2*rim+1)))
bkgAB = np.empty(((OB.nfiles/2,2*rbg+1,2*rbg+1)))
fluxAB = np.empty(OB.nfiles/2)
for i in range(OB.nfiles/2):
    iA, iB = 2*i, 2*i+1
    cx, cy = OB.cx[iA], OB.cy[iA]
    
    # first case: the off-axis profile, 2 images per acquisition (chop A and chop B)
    if targetNum == 0:
        
        # open fits and load data
        psfA = OB.getData(iA, 1) # chop A
        psfB = OB.getData(iB, 1) # chop B
        diff = (psfB-psfA)/2. # chopped image
        
        # integrate the psf signal within the FWHM (circular aperture photometry)
        psfAB[i,:,:] = diff[cy-rim:cy+rim+1,cx-rim:cx+rim+1]
        r = ImPro.get_r_dist(2*rim+1,2*rim+1,rim,rim)
        circ_area = np.where(r < OB.FWHM[iA]/2.)
        fluxAB[i] = np.sum(psfAB[i,:,:][circ_area])
        #plt.figure(); plt.imshow(psfAB[i,:,:], interpolation='nearest', cmap="CMRmap"); plt.colorbar()
        
        # fit an Airy pattern to check the off-axis PSF centers
        if offaxisFit == True:
            try:
                (A,xo,yo,F) = ImPro.fit_airy_2D(psfAB[i,:,:])
            except:
                xo,yo = 0,0
            OB.cx2.append(cx-rim+xo)
            OB.cy2.append(cy-rim+xo)
    
    # other cases: on sky, 2 images per acquisition (nod A and nob B, pre-chopped)
    else:
        
        # open fits and load data ([3] is the chopped image)
        imgA = OB.getData(iA, 3) # nod A
        imgB = OB.getData(iB, 3) # nod B
        
        # fit a gaussian to find the exact centers
        xstart, xend = cx-rim, cx+rim+1
        ystart, yend = cy-rim, cy+rim+1
        psfA = imgA[ystart-chopnod:yend-chopnod, xstart:xend]
        psfB = imgB[ystart+chopnod:yend+chopnod, xstart:xend]
        (A,xa,ya,F) = ImPro.fit_airy_2D(-psfA)
        (A,xb,yb,F) = ImPro.fit_airy_2D(psfB)
        OB.cx2.append(cx + (xa+xb-2*rim)/2.)
        OB.cy2.append(cy + (ya+yb-2*rim)/2.)
        
        # build the attenuated PSF (AGPM on-axis)
        cx2 = np.round(OB.cx2[-1]).astype(int)
        cy2 = np.round(OB.cy2[-1]).astype(int)
        psfA_AGPM = imgA[cy2-rim:cy2+rim+1,cx2-rim:cx2+rim+1]
        psfB_AGPM = imgB[cy2-rim:cy2+rim+1,cx2-rim:cx2+rim+1]
        psfAB_AGPM[i,:,:] = (psfB_AGPM-psfA_AGPM)/2.
        
        # for chop-nod mode, build the off-axis PSF, and the background
        chopnod2 = chopnod + (yb-ya)/2
        cy2A = np.round(OB.cy2[-1]-chopnod2).astype(int)
        cy2B = np.round(OB.cy2[-1]+chopnod2).astype(int)
        psfA = imgA[cy2A-rim:cy2A+rim+1,cx2-rim:cx2+rim+1]
        psfB = imgB[cy2B-rim:cy2B+rim+1,cx2-rim:cx2+rim+1]
        psfAB[i,:,:] = (psfB-psfA)/2.
        bkgA = imgA[cy2-rbg:cy2+rbg+1,cx2-rbg:cx2+rbg+1]
        bkgB = imgB[cy2-rbg:cy2+rbg+1,cx2-rbg:cx2+rbg+1]
        bkgAB[i,:,:] = (bkgB-bkgA)/2.
        
        #psfAB_sub = psfAB - np.median(bkgAB)   #-np.median(psfAB)
        #psfAB_AGPM_sub = psfAB_AGPM + np.median(bkgAB)   #-np.median(psfAB_AGPM)




#######



""" Figures """
""" ******* """ 

# 1/ Off-axis Transmission

if targetNum == 0:
    sequence = np.array([13,11,10,9,8,7,6,2,1,0]);
    fsz = 7#9; #fonsize
    lim = rim/OB.lamD[0]
    #rim2 = int(2*OB.FWHM[0]);
    
    # LINEAR SEQUENCE
    vmin = 0.; vmax = 1.;
    f=plt.figure();
#    f.subplots_adjust(wspace=0)
    f.subplots_adjust(wspace=0.,right=.8);
    cbar_ax = f.add_axes([.81, 0.455, 0.007, 0.088]);
    # offsets in lamD
    off = pxmm/OB.lamD[0]*(114.50-np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05]))
    for k,i in enumerate(sequence): 
        myAxes = [-lim+off[i],lim+off[i],-lim,lim];
#        psf = PSFs[i,cim-rim2:cim+rim2+1,cim+offsets[i]-rim2:cim+offsets[i]+rim2+1]/np.max(PSFs);
        psf = psfAB[i]/np.max(psfAB)
        ax = plt.subplot2grid((1,sequence.size),(0,k));
        im = ax.imshow(psf,origin='lower',vmin=vmin,vmax=vmax,extent=myAxes,\
                interpolation='nearest',cmap="CMRmap");
        ax.tick_params(labelsize=fsz,direction='out',length=3,right='off',top='off');
        ax.tick_params(left='off');
        ax.set_xticks([round(off[i],1)]);
        ax.set_yticks([-1,0,1]);  # is bugged, should appear on the figure...
        if k == 0:
            ax.tick_params(left='on');
            ax.set_ylabel(r"$\lambda$/D",fontsize=fsz)
            cb = f.colorbar(im, cax=cbar_ax);
            cb.ax.tick_params(labelsize=fsz);
            cb.set_ticks([vmin,(vmax+vmin)/2.,vmax]);
        if k == round(sequence.size/2):
            ax.set_xlabel(r"$\lambda$/D",fontsize=fsz);
    plt.setp([a.get_yticklabels() for a in f.axes[1:]], visible=False);   
    plt.savefig(os.path.join(path,'offaxislin.png'), dpi=300, transparent=True)
    
    # LOG10 SEQUENCE
    vmin = -4;vmax = 0;
    f=plt.figure();
    f.subplots_adjust(wspace=0.,right=.8);
    cbar_ax = f.add_axes([.81, 0.455, 0.007, 0.088]);
    # offsets in lamD
    off = pxmm/OB.lamD[0]*(114.50-np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05]))
    for k,i in enumerate(sequence): 
        myAxes = [-lim+off[i],lim+off[i],-lim,lim];
#        psf = PSFs[i,cim-rim2:cim+rim2+1,cim+offsets[i]-rim2:cim+offsets[i]+rim2+1]/np.max(PSFs);
        psf = psfAB[i]/np.max(psfAB)
        psf = np.log10(psf); # for log10
        psf[np.where(np.isnan(psf))] = -4 # set nan values to 10^-4
        ax = plt.subplot2grid((1,sequence.size),(0,k));
        im = ax.imshow(psf,origin='lower',vmin=vmin,vmax=vmax,extent=myAxes,\
                interpolation='nearest',cmap="CMRmap");
        ax.tick_params(labelsize=fsz,direction='out',length=3,right='off',top='off');
        ax.tick_params(left='off');
        ax.set_xticks([round(off[i],1)]);
        ax.set_yticks([-1,0,1]);  # is bugged, should appear on the figure...
        ax.add_artist(plt.Circle((off[i],0),0.5,color='w',fill=False));
        if k == 0:
            ax.tick_params(left='on');
            ax.set_ylabel(r"$\lambda$/D",fontsize=fsz)
            cb = f.colorbar(im, cax=cbar_ax);
            cb.ax.tick_params(labelsize=fsz);
            cb.set_ticks([vmin,(vmax+vmin)/2.,vmax]);
        if k == round(sequence.size/2):
            ax.set_xlabel(r"$\lambda$/D",fontsize=fsz);
    plt.text(-1,2.3,"log10",fontsize=fsz)
    plt.text(-39.,2.3,"circles = FWHM",fontsize=fsz)
    plt.setp([a.get_yticklabels() for a in f.axes[1:]], visible=False);   
    plt.savefig(os.path.join(path,'offaxislog.png'), dpi=300, transparent=True)
    
    # RADIAL PROFILE
    f=plt.figure();
    plt.xlabel(r"$\lambda$/D")
    plt.ylabel("Radial profile")
    x1=pxmm/OB.lamD[0]*(np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05])-114.50)
    y1=fluxAB/np.max(fluxAB)*0.838 # 0.838 = facteur d'Elsa
    x2 = np.array(sorted(abs(x1)))
    y2 = np.array(sorted(y1))
    dat = np.loadtxt(os.path.join(instruPath, 'Elsa/VLT_circEP_Lyot_circ_aperture_nobsc_lyotout90_lyotin20_lbd37_custom-norm.txt'))
    x3 = dat[:,0]
    y3 = dat[:,1]
    plt.grid('on')
#    plot([-3,3], [.5,.5], 'k--')
    plt.xlim([-3,3])
    plt.plot(x3,y3, 'g', label="simulated")
    plt.plot(x2,y2, '+m', markersize=8, mew=2, label="measured")
    plt.yscale('log')
    plt.plot([0,3], [.5,.5], 'k--')
    plt.xlim([0,3])
    plt.ylim(2e-2,1.5)
    IWA = 1.14 # 1.14 de Elsa
    plt.plot([IWA,IWA], [2e-2,.5], 'k--', label='IWA = '+str(round(IWA,2))+'$\lambda$/D')
    plt.xlabel(r"$\lambda$/D")
    plt.legend(loc='lower right', fontsize=8)
    plt.savefig(os.path.join(path,'radprof.png'), dpi=300, transparent=True)
    
    # xo,yo position fit
    if offaxisFit == True:
        dist = np.array([offsets[2*i] for i in range(OB.nfiles/2)])
        x = np.array([OB.cx[2*i] for i in range(OB.nfiles/2)])
        y = np.array([OB.cy[2*i] for i in range(OB.nfiles/2)])
        ymed = np.median(OB.cy2)
        # xo position
        plt.figure();
        plt.plot([dist[0],dist[-1]],[x[0],x[-1]])
        plt.plot(dist,x,'g+')
        plt.plot(dist[:3],OB.cx2[:3], 'r+-', label="Off-axis X-distance")
        plt.plot(dist[-7:],OB.cx2[-7:], 'r+-', label="Off-axis X-distance")
        plt.xlim(114,115); plt.ylim(370,410)
        plt.xlabel("x offset (mm)"); plt.ylabel("x position (pixels)")
        plt.savefig(os.path.join(path,'xo-position.png'), dpi=300); plt.clf()
        # yo position
        plt.figure();
        plt.plot([dist[0],dist[-1]],[y[0],y[0]])
        plt.plot(dist,y,'g+')
        plt.plot(dist[:3],OB.cy2[:3], 'r+-', label="Off-axis Y-distance")
        plt.plot(dist[-6:],OB.cy2[-6:], 'r+-', label="Off-axis Y-distance")
        plt.xlim(114,115); plt.ylim(160,170)
        plt.xlabel("x offset (mm)"); plt.ylabel("y position (pixels)")
        plt.savefig(os.path.join(path,'yo-position.png'), dpi=300); plt.clf()
# OB.cx, OB.cy, OB.cx2, OB.cy2 


# 2/ On-sky targets

#   myAxes = np.array([-1,1,-1,1])*rim/OB.lamD[0]
else:
    f=plt.figure();
    f.subplots_adjust(hspace=0)
#    f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=14)
    ax1 = plt.subplot2grid((2,2), (0, 0))
    ax2 = plt.subplot2grid((2,2), (1, 0))
    ax3 = plt.subplot2grid((2,2), (0, 1),rowspan=2)
    ax1.imshow(imgB[cy-rbg/3.:cy+rbg+1,cx-rbg*.7:cx+rbg*.7+1],origin='lower',vmin=np.min(bkgAB),vmax=np.max(bkgAB));
    ax2.imshow(-imgA[cy-rbg:cy+rbg/3+1.,cx-rbg*.7:cx+rbg*.7+1],origin='lower',vmin=np.min(bkgAB),vmax=np.max(bkgAB)); 
    im=ax3.imshow(bkgAB,origin='lower');
    ax3.set_title('Nod B - Nod A \n Bkg median='+str(round(np.median(bkgAB),2)))
    plt.setp([a.get_xticklabels() for a in f.axes[:-2]], visible=False)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    f.subplots_adjust(right=0.8)
    f.colorbar(im, cax=cbar_ax)
    plt.savefig(path+filt+'_bkg.png', dpi=300, transparent=True);plt.clf()
    plt.savefig(os.path.join(path,'offaxislog.png'), dpi=300, transparent=True)



    
    Flux_nocor_circ=  np.sum(psfAB[circ_area])
    Flux_AGPM_circ = np.sum(psfAB_AGPM[circ_area]) 
    N = round(Flux_nocor_circ/Flux_AGPM_circ,1)
    print '\nAperture photometry:'
    print('Rejection = '+str(N))
    " Radial profiles "
    P_nocor = ImPro.get_radial_profile(psfAB, (xo,yo), nbin, disp=0)
    P_AGPM = ImPro.get_radial_profile(psfAB_AGPM, (xo,yo), nbin, disp=0)
    # normalization
    norm = np.max(P_nocor) 
    P_nocor = P_nocor/norm
    P_AGPM  = P_AGPM/norm
    # constrained minimal value
    P_nocor = P_nocor - np.minimum(0,np.min(P_nocor)-minprofiles)
    P_AGPM = P_AGPM - np.minimum(0,np.min(P_AGPM)-minprofiles)
    # corresponding radial distance
    (nxp) = P_nocor.shape[0]
    Xlbd=np.linspace(0,nxp-1,nxp)/lamD 

    f=plt.figure();
    f.subplots_adjust(hspace=0)
    f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=14)
    ax1 = plt.subplot2grid((2,4), (0, 0),rowspan=2,colspan=3)
    ax2 = plt.subplot2grid((2,4), (0, 3))
    ax3 = plt.subplot2grid((2,4), (1, 3))
    ax1.plot(Xlbd*nbin, P_nocor, 'r', label="Off-axis PSF")
    ax1.plot(Xlbd*nbin, P_AGPM, 'b', label="Attenuated PSF \n R="+str(N)+" ; seeing="+str(seeing)+"''")
    ax1.set_ylabel("Normalized PSF radial profile")
    ax1.set_xlabel(r"$\lambda$/D")
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_yscale('log'); ax1.set_ylim(minprofiles,1.);ax1.set_xlim(0,5.)
#        ax1.set_title('Radial profiles')
    im=ax2.imshow(psfAB,origin='lower',extent=myAxes); 
    ax3.imshow(psfAB_AGPM,origin='lower',extent=myAxes,vmin=np.min(psfAB),vmax=np.max(psfAB)); 
    ax3.set_xlabel(r"$\lambda$/D")
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    f.subplots_adjust(right=0.8)
    f.subplots_adjust(top=0.9)
    f.colorbar(im, cax=cbar_ax)
    # markers (FWHM, PSF zeros)
    R=OB.FWHM[iA]/OB.lamD[iA]/2.; ax1.plot([R,R], [minprofiles, 1.], 'k--')
    #R=1*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
    #R=2*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
    #R=3*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
    plt.savefig(path+filt+'_radprof.png', dpi=300, transparent=True);plt.clf()



