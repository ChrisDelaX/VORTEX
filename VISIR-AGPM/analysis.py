# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 07:11:43 2015

@author: cdelacroix
"""

import os.path
import glob as glob
import numpy as np
import PSFtools as PSF

import matplotlib.pyplot as plt
import scipy as sp
from astropy.io import fits         # fits file reading/writing


""" Input Data """

instruPath = os.path.expandvars('$HOME/INSTRUMENTS/VISIR')
minprofiles = 1e-4  # constrained minimal value for the radial values, to allow log scale
nbin = 1.           # binning parameter for the radial profile
asecpmm = 1.9484    # cf. Eric Pantin
asecpp = 45.3e-3    # cf. Eric Pantin (arcsec per pixel = asecpmm*.8/34)
LS = 0.9            # Lyot-stop
diam = 8.2          # telescope diameter
chopnod = 178       # chopping-nodding offset (in pixels)
rim = 50            # image region of interest
rbg = 250           # background region of interest
targetNum = 0       # select the set of images corresponding to a run/target
offaxisFit = True   # to find the off-axis PSF centers


""" Targets and filters """
offsets = []
if targetNum == 0:
    target = '2015-03-09_off-axis'
    images = ['0020','0021','0023','0022','0026','0027','0030','0031','0033','0032','0034','0035','0037','0036','0038','0039','0041','0040','0042','0043','0044','0045','0046','0047','0048','0049','0050','0052','0053','0054']
    offsets = [114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05] # from Eric
    # to print offsets as string: ('%.2f' % offsets[i]).zfill(6)
    nfits = np.size(offsets)
    filters = ['BB-AGPM']*nfits
    lams = np.array([12.4]*nfits)*1e-6
    seeings = [0]*nfits
    rim = 15 # specific to this run
    chopnod = 0 # specific to this run
elif targetNum == 1:
    target = 'HD125932'
    images = ['0009','0010','0011','0012','0013','0014','0015','0016','0017','0018','0019','0020']
    filters = ['BB-AGPM','PAH1','NB-4QPM1','NeII','PAH2','NB-4QPM2']
    lams = np.array([12.4,8.6,10.65,12.4,11.3,11.3])*1e-6
    seeings=([1.07,1.21,1.16,1.70,1.18,0.84])
elif targetNum == 2:
    target = 'HD139063'
    images = ['0021','0022','0023','0024','0025','0026','0027','0028','0029','0030']
    filters = ['BB-AGPM','PAH1','NB-4QPM1','NB-4QPM2','BB-AGPMbis']
    lams = np.array([12.4,8.6,10.65,11.3,12.4])*1e-6
    seeings=([0.89,0.90,0.72,0.83,0.96])
    bws = np.array([1.6,1.,0,0,1.6])*1e-6 # ????
elif targetNum == 3:
    target = 'HR6378'
    images = ['0031','0032']
    filters = ['BB-AGPM']
    lams = np.array([12.4])*1e-6
    seeings=([1.22])

path = os.path.join(instruPath, target)
nfits = np.size(filters)
offsets = [0]*nfits if not len(offsets) else offsets
ofaxMat = np.empty(((nfits,2*rim+1,2*rim+1)))
fluxMat = np.empty(nfits)


""" Off-axis profile """

if targetNum==0:
    
    # find the off-axis PSF centers
    if offaxisFit == True:
        dist = [114.05,114.15,114.2,114.25,114.3,114.35,114.4,114.65,114.75,114.85]
        x = [373.39,377.54,379.64,381.57,383.32,384.72,384.46,400.18,403.58,407.72] # from Gaussian fit
        y = [164.43,164.45,164.53,164.55,164.61,164.58,164.67,164.3,163.72,163.51] # from Gaussian fit
        ymed = np.median(y[:7])
        dist2 = np.zeros(np.size(offsets))
        x2 = np.zeros(np.size(offsets))
        y2 = np.zeros(np.size(offsets))
        for k in range(np.size(offsets)):
            dist2[k] = offsets[k]
            x2[k] = 372.89+.429125*(1+100*(offsets[k]-114.05))
            y2[k] = 164.55 #ymed    

        plt.figure(1); plt.clf()
        plt.plot([dist[0],dist[-1]],[x[0],x[-1]])
        plt.plot(dist[:7],x[:7], 'r+-', label="Off-axis X-distance")
        plt.plot(dist[7:],x[7:], 'r+-', label="Off-axis X-distance")
        plt.plot(dist2,x2,'g+')
        plt.xlabel("x offset (mm)"); plt.ylabel("x position (pixels)")
        plt.savefig('xo-position.png', dpi=300);plt.clf()


        plt.figure(2); plt.clf()
        #plt.plot([dist[0],dist[-1]],[y[0],y[0]])
        plt.plot([dist[0],dist[-1]],[y2[0],y2[0]])
        plt.plot(dist[:7],y[:7], 'r+-', label="Off-axis Y-distance")
        plt.plot(dist[7:],y[7:], 'r+-', label="Off-axis Y-distance")
        plt.ylim(160,170)
        plt.plot(dist2,y2,'g+')
        plt.xlabel("x offset (mm)"); plt.ylabel("y position (pixels)")
        plt.savefig('yo-position.png', dpi=300);plt.clf()






rim=50


cx,cy = [391, 437]     # obs 09-10,13-14
#for i,(filt,lb,seeing) in enumerate(zip(filters,lams,seeings)):

for i in range(1):
    obsA=images[2*i]
    obsB=images[2*i+1]
    if obsA in {'0011','0015','0017','0023'} or obsB in {'0012','0016','0018','0024'}:
        cy = 549 # obs 11-12,15-18, 23-24
    if targetNum == 0:
        cx = int(np.round(372.89+42.9125*(offsets[i]-114.05))) # cf. fit off-axis distance
        cy = 164
    disp = 1
    
    # open fits and load data
    pathA = glob.glob(os.path.join(path, '*'+obsA+'*.fits'))[0]
    pathB = glob.glob(os.path.join(path, '*'+obsB+'*.fits'))[0]
    if chopnod == 0: # off-axis transmission sequence
        with fits.open(pathA) as hdulist:
            header0 = hdulist[0].header
            header1 = hdulist[1].header
            filename = header0.get('ORIGFILE')
            assert pathA == os.path.join(path, filename), 'wrong filename'
            # save header
            try:
                with open(filename[:-4]+'txt', 'r') as f:
                    pass
            except IOError:
                with open(filename[:-4]+'txt', 'w') as f:
                    header = open(filename, 'w')
            
            psfA = hdulist[1].data
            
            
            
        with fits.open(pathB) as hdulist:
            header0 = hdulist[0].header
            header1 = hdulist[1].header
            filename = header0.get('ORIGFILE')
            assert pathB == os.path.join(path, filename), 'wrong filename'
            psfB = hdulist[1].data
            # filename = hdulist[0].header.get('ORIGFILE')
            # expTime = hdulist[1].header.get('EXPTIME')    # 0.0024887
        PSF_AGPM = (PSF_AGPM_B-PSF_AGPM_A)/2.


import ImageProcessing, ImageProcessing.ObservingBlock

reload(ImageProcessing);reload(ImageProcessing.ObservingBlock)
from ImageProcessing.ObservingBlock import ObservingBlock
folder = '/Users/cdelacroix/INSTRUMENTS/VISIR/2015-03-09_off-axis'
OB = ObservingBlock(folder, start='corono')
print OB



# HIERARCH ESO DET CHOP TIM    = 0.0024887 / [s] TIM period
# HIERARCH ESO DET CHOP TRANSTIM= 0.025 / [s] Chopper transition time
# HIERARCH ESO DET SEQ1 DIT    = 0.0024887 / [s] Integration time
# HIERARCH ESO DET SEQ1 INITTIME= 0.0000000 / [s] Exposure Init Time
# HIERARCH ESO DET SEQ1 EXPTIME= 12.4460887 / [s] Exposure Sequence Time
# HIERARCH ESO DET SEQ1 MINDIT = 0.0024887 / [s] Minimum DIT
# 
        
    else: # target observations: nod A and nod B (already chopped)
        with fits.open(pathA) as hdulist:
            psfA = hdulist[3].data # [3] is the chopped image
        with fits.open(pathB) as hdulist:
            psfB = hdulist[3].data # [3] is the chopped image
        
        # fit a gaussian to find the exact centers
        PSF_nocor_A = psfA[cy-rim-chopnod:cy+rim-chopnod+1,cx-rim:cx+rim+1]
        PSF_nocor_B = psfB[cy-rim+chopnod:cy+rim+chopnod+1,cx-rim:cx+rim+1]
        (A,xa,ya,F) = PSF.fit_airy_2D(-PSF_nocor_A)
        (A,xb,yb,F) = PSF.fit_airy_2D(PSF_nocor_B)
        cx2 = cx + (xa+xb-2*rim)/2
        cy2 = cy + (ya+yb-2*rim)/2
        chopnod2 = chopnod + (yb-ya)/2
    
    
    # build the attenuated PSF (AGPM on-axis)
    PSF_AGPM_A = psfA[cy2-rim:cy2+rim+1,cx2-rim:cx2+rim+1]
    PSF_AGPM_B = psfB[cy2-rim:cy2+rim+1,cx2-rim:cx2+rim+1]
    PSF_AGPM = (PSF_AGPM_B-PSF_AGPM_A)/2.
    
    # for chop-nod mode, build the off-axis PSF, and the background
    PSF_nocor_A = psfA[cy2-rim-chopnod2:cy2+rim-chopnod2+1,cx2-rim:cx2+rim+1]
    PSF_nocor_B = psfB[cy2-rim+chopnod2:cy2+rim+chopnod2+1,cx2-rim:cx2+rim+1]
    PSF_nocor = (PSF_nocor_B-PSF_nocor_A)/2.
    bkg_A = psfA[cy2-rbg:cy2+rbg+1,cx2-rbg:cx2+rbg+1]
    bkg_B = psfB[cy2-rbg:cy2+rbg+1,cx2-rbg:cx2+rbg+1]
    bkg = (bkg_B-bkg_A)/2.

    plt.figure()
    plt.imshow(PSF_AGPM_A)
    colorbar()
    plt.figure()
    plt.imshow(PSF_AGPM_B)
    colorbar()

    plt.figure()
    plt.imshow(PSF_AGPM_B-PSF_AGPM_A)
    colorbar()
    
    try:
        (A,xb,yb,F) = PSF.fit_airy_2D(PSF_AGPM)
    except:
        xa,ya, xb,yb = 0,0,0,0

    print xa,xb,ya,yb


    #PSF.ensure_dir(path+filt+'/')


    PSF_nocor[PSF_nocor<minprofiles] = minprofiles   
    PSF_AGPM[PSF_AGPM<minprofiles] = minprofiles   
    #PSF_nocor_sub = PSF_nocor - np.median(bkg)   #-np.median(PSF_nocor)
    #PSF_AGPM_sub = PSF_AGPM + np.median(bkg)   #-np.median(PSF_AGPM)
    
    " Circular aperture photometry " 
    lamD = lb/(diam*LS)*(60*60*180/np.pi) # lam/D in arcsec
    lbdDpp = asecpp/lamD  # lb/d per pixel
    FWHM = 1/lbdDpp # Full Width at Half Maximum = lambda/D
    myAxes = [-rim*lbdDpp,rim*lbdDpp,-rim*lbdDpp,rim*lbdDpp]
    (xo,yo) = (rim,rim)
    (nx,ny)=(2*rim+1,2*rim+1)
    Rcirc = FWHM/2.           # radius of the integrated zone
    r = PSF.get_r_dist(nx,ny,xo,yo)
    circ_area = np.where(r<Rcirc)
    
    if targetNum == 0:
        ofaxMat[i,:,:] = PSF_nocor
        fluxMat[i] = np.sum(PSF_nocor[circ_area])
        if i == (nfits-1):
            fluxNorm = fluxMat[i] 
            
    else:
        f=plt.figure(1);
        f.subplots_adjust(hspace=0)
        f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=14)
        ax1 = plt.subplot2grid((2,2), (0, 0))
        ax2 = plt.subplot2grid((2,2), (1, 0))
        ax3 = plt.subplot2grid((2,2), (0, 1),rowspan=2)
        ax1.imshow(psfB[cy-rbg/3.:cy+rbg+1,cx-rbg*.7:cx+rbg*.7+1],origin='lower',vmin=np.min(bkg),vmax=np.max(bkg));
        ax2.imshow(-psfA[cy-rbg:cy+rbg/3+1.,cx-rbg*.7:cx+rbg*.7+1],origin='lower',vmin=np.min(bkg),vmax=np.max(bkg)); 
        im=ax3.imshow(bkg,origin='lower');
        ax3.set_title('Nod B - Nod A \n Bkg median='+str(round(np.median(bkg),2)))
        plt.setp([a.get_xticklabels() for a in f.axes[:-2]], visible=False)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.subplots_adjust(right=0.8)
        f.colorbar(im, cax=cbar_ax)
        plt.savefig(path+filt+'_bkg.png', dpi=300, transparent=True);plt.clf()
    
        Flux_nocor_circ=  np.sum(PSF_nocor[circ_area])
        Flux_AGPM_circ = np.sum(PSF_AGPM[circ_area]) 
        N = round(Flux_nocor_circ/Flux_AGPM_circ,1)
        print '\nAperture photometry:'
        print('Rejection = '+str(N))
        " Radial profiles "
        P_nocor = PSF.get_radial_profile(PSF_nocor, (xo,yo), nbin, disp=0)
        P_AGPM = PSF.get_radial_profile(PSF_AGPM, (xo,yo), nbin, disp=0)
        # normalization
        norm = np.max(P_nocor) 
        P_nocor = P_nocor/norm
        P_AGPM  = P_AGPM/norm
        # constrained minimal value
        P_nocor = P_nocor - np.minimum(0,np.min(P_nocor)-minprofiles)
        P_AGPM = P_AGPM - np.minimum(0,np.min(P_AGPM)-minprofiles)
        # corresponding radial distance
        (nxp) = P_nocor.shape[0]
        Xlbd=np.linspace(0,nxp-1,nxp)*lbdDpp 
    
        f=plt.figure(2);
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
        im=ax2.imshow(PSF_nocor,origin='lower',extent=myAxes); 
        ax3.imshow(PSF_AGPM,origin='lower',extent=myAxes,vmin=np.min(PSF_nocor),vmax=np.max(PSF_nocor)); 
        ax3.set_xlabel(r"$\lambda$/D")
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.subplots_adjust(right=0.8)
        f.subplots_adjust(top=0.9)
        f.colorbar(im, cax=cbar_ax)
        # markers (HWFM, PSF zeros)
        R=Rcirc*lbdDpp; ax1.plot([R,R], [minprofiles, 1.], 'k--')
        #R=1*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
        #R=2*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
        #R=3*1.22; ax1.plot([R,R], [minprofiles, 1.], 'm:')
        plt.savefig(path+filt+'_radprof.png', dpi=300, transparent=True);plt.clf()

            
""" Off-axis Transmission """
""" ********************* """ 

if targetNum == 0:
    f=plt.figure(1);
    f.subplots_adjust(wspace=0)
    #f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=8)
    ax1 = plt.subplot2grid((1,12), (0, 0));im1=ax1.imshow(ofaxMat[14],origin='lower',extent=myAxes); 
    ax1.set_ylabel(r"$\lambda$/D",fontsize=7)
    ax2 = plt.subplot2grid((1,12), (0, 1));ax2.imshow(ofaxMat[13],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax3 = plt.subplot2grid((1,12), (0, 2));ax3.imshow(ofaxMat[12],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax4 = plt.subplot2grid((1,12), (0, 3));ax4.imshow(ofaxMat[11],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax5 = plt.subplot2grid((1,12), (0, 4));ax5.imshow(ofaxMat[10],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax6 = plt.subplot2grid((1,12), (0, 5));ax6.imshow(ofaxMat[9],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax7 = plt.subplot2grid((1,12), (0, 6));ax7.imshow(ofaxMat[8],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax7.set_xlabel(r"$\lambda$/D",fontsize=7)
    ax8 = plt.subplot2grid((1,12), (0, 7));ax8.imshow(ofaxMat[7],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax9 = plt.subplot2grid((1,12), (0, 8));ax9.imshow(ofaxMat[6],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax13 = plt.subplot2grid((1,12), (0, 9));ax13.imshow(ofaxMat[2],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax14 = plt.subplot2grid((1,12), (0, 10));ax14.imshow(ofaxMat[1],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    ax15 = plt.subplot2grid((1,12), (0, 11));ax15.imshow(ofaxMat[0],origin='lower',extent=myAxes,vmin=np.min(ofaxMat[0]),vmax=np.max(ofaxMat[0])); 
    plt.setp([a.get_yticklabels() for a in f.axes[1:]], visible=False)
    for a in f.axes:
        a.tick_params(labelsize=5,direction='out',length=3)
        a.set_yticks([-1,0,1])
        a.set_xticks([-1,0,1])
    cbar_ax = f.add_axes([0.82, 0.45, 0.01, 0.1])
    f.subplots_adjust(right=0.8)
    cb = f.colorbar(im1, cax=cbar_ax)
    cb.set_ticks([0,20,40])
    plt.savefig(path+filt+'_PSFs.png', dpi=300, transparent=True);#plt.clf()
    
    
    f=plt.figure(2);
    x1=asecpmm/lamD*(np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05])-114.50)
    y1=fluxMat/fluxNorm
    ax1 = plt.subplot2grid((1,2), (0, 0));
    ax1.plot(x1,fluxMat/fluxNorm, 'r+--', label="Off-axis PSF")
    ax1.grid('on')
    ax1.plot([-3,3], [.5,.5], 'k--')
    ax1.set_xlim([-3,3])
    ax1.set_xlabel(r"$\lambda$/D")
    ax1.set_ylabel("Off-axis transmission")
    x2 = np.array(sorted(abs(x1)))
    y2 = np.array(sorted(y1))
    x3 = np.linspace(.1,3.,25)
    y3 = 1 - 4*sp.special.j1(np.pi*x3/np.sqrt(3))**2 / (np.pi*x3/np.sqrt(3))**2
    [a,b,c] = PSF.fit_gauss_1D(-(y2-1),x2)
    y4 = -(PSF.oneD_Gaussian(x3, a, b, c)-1)
    ax2 = plt.subplot2grid((1,2), (0, 1));
    ax2.plot(x2,y2, 'r+', markersize=4, label="measures")
    ax2.plot(x3,y3,'g', label="Riaud&Hanot10")
    ax2.plot(x3,y4,'k', label="Fitted Gaussian")
    ax2.set_yscale('log')
    ax2.grid('on')
    ax2.set_xlim([0,3])
    ax2.set_ylim(2e-2,1.5)
    ax2.plot([0,3], [.5,.5], 'k--')
    
    IWA = np.sqrt(2*c**2*np.log(2*a)) + b
    ax2.plot([IWA,IWA], [2e-2,.5], 'b-.', label='IWA = '+str(round(IWA,2))+'$\lambda$/D')
    ax2.set_xlabel(r"$\lambda$/D")
    ax2.legend(loc='lower right', fontsize=8)
    
    plt.savefig(path+filt+'_offaxis.png', dpi=300, transparent=True);#plt.clf()


    
    




