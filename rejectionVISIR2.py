# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 07:11:43 2015

@author: cdelacroix
"""

import scipy as sp
import numpy as np
import AGPM_lib_visir as AGPM
import matplotlib.pyplot as plt     


""" Input Data """
""" ********** """

instru = 'VISIR'
minprofiles = 1e-4 # constrained minimal value for the radial values, to allow log scale
nbin = 1. # binning parameter for the radial profile

obsSelect = 0
if obsSelect==0:
    target = '2015-03-09_off-axis'
    obsnum = ['0020','0021','0023','0022','0026','0027','0030','0031','0033','0032','0034','0035','0037','0036','0038','0039','0041','0040','0042','0043','0044','0045','0046','0047','0048','0049','0050','0052','0053','0054']
    xOffset = ['114.85','114.75','114.65','114.55','114.50','114.50','114.50','114.45','114.40','114.35','114.30','114.25','114.20','114.15','114.05']    
    nfits = np.size(xOffset)
    filters = ['BB-AGPM']*nfits
    lbs = np.array([12.4]*nfits)*1e-6
    seeings = [0]*nfits
    offset=0    
    rim = 15
elif obsSelect==1:
    target = 'HD125932'
    obsnum = ['0009','0010','0011','0012','0013','0014','0015','0016','0017','0018','0019','0020']
    filters = ['BB-AGPM','PAH1','NB-4QPM1','NeII','PAH2','NB-4QPM2']
    nfits = np.size(filters)
    lbs = np.array([12.4,8.6,10.65,12.4,11.3,11.3])*1e-6
    seeings=([1.07,1.21,1.16,1.70,1.18,0.84])
    xOffset = [0]*nfits
    offset = -178.
    rim = 50.
elif obsSelect==2:
    target = 'HD139063'
    obsnum = ['0021','0022','0023','0024','0025','0026','0027','0028','0029','0030']
    filters = ['BB-AGPM','PAH1','NB-4QPM1','NB-4QPM2','BB-AGPMbis']
    nfits = np.size(filters)
    lbs = np.array([12.4,8.6,10.65,11.3,12.4])*1e-6
    bws = np.array([1.6,1.,0,0,1.6])*1e-6
    seeings=([0.89,0.90,0.72,0.83,0.96])
    xOffset = [0]*nfits
    offset = -178.
    rim = 50.
elif obsSelect==3:
    target = 'HR6378'
    obsnum = ['0031','0032']
    filters = ['BB-AGPM']
    nfits = np.size(filters)  
    lbs = np.array([12.4])*1e-6
    seeings=([1.22])
    xOffset = [0]*nfits
    offset = -178.
    rim = 50.
path = '/Users/cdelacroix/INSTRUMENTS/'+instru+'/'+target+'/'


""" Region of interest """
rbkg = int(250)
asecpmm = 1.9484  #says Eric
asecpp = 45.3e-3  #says Eric   # arcsec per pixel =asecpmm*.8/34
LS = 0.9
diam = 8.2

ofaxMat = np.empty(((nfits,2*rim+1,2*rim+1)))
fluxMat = np.empty(nfits)
cxx = np.empty(rim)
for k in range(nfits):
    obsA=obsnum[2*k]
    obsB=obsnum[2*k+1]
    cx=391.;cy=437.; # obs 09-10,13-14
    if (obsA=='0011') or (obsB=='0012') or (obsA=='0015') or (obsB=='0016') or (obsA=='0017') or (obsB=='0018') or (obsA=='0023') or (obsB=='0024'):
        cy=549.; # obs 11-12,15-18, 23-24
    if obsSelect == 0:
        cx= 372.89+.429125*(int(xOffset[k][-2:])-5)
        cy= 164.
    center=[int(round(cx)),int(cy)]
    filt=filters[k]
    lb=lbs[k]
    seeing=seeings[k]
    disp=1
    #AGPM.ensure_dir(path+filt+'/')
    (bkg,PSF_nocor,PSF_AGPM,HDU_A,HDU_B) = AGPM.get_image_VISIR(k,path,obsA,obsB,center,rim,rbkg,filt,target,lb,seeing,offset,xOffset[k])
    PSF_nocor[PSF_nocor<minprofiles] = minprofiles   
    PSF_AGPM[PSF_AGPM<minprofiles] = minprofiles   
    #PSF_nocor_sub = PSF_nocor - np.median(bkg)   #-np.median(PSF_nocor)
    #PSF_AGPM_sub = PSF_AGPM + np.median(bkg)   #-np.median(PSF_AGPM)
    
    " Circular aperture photometry " 
    lb = lbs[k]
    asecplbd = lb/(diam*LS)*(60*60*180/np.pi) #arcsec per lbD
    lbdDpp = asecpp/asecplbd  # lb/d per pixel
    FWHM = 1/lbdDpp # Full Width at Half Maximum = lambda/D
    myAxes = [-rim*lbdDpp,rim*lbdDpp,-rim*lbdDpp,rim*lbdDpp]
    (xo,yo) = (rim,rim)
    (nx,ny)=(2*rim+1,2*rim+1)
    Rcirc = FWHM/2.           # radius of the integrated zone
    r = AGPM.get_r_dist(nx,ny,xo,yo)
    circ_area = np.where(r<Rcirc)
    
    if obsSelect == 0:
        ofaxMat[k,:,:] = PSF_nocor
        fluxMat[k] = np.sum(PSF_nocor[circ_area])
        if k ==(nfits-1):
            fluxNorm = fluxMat[k] 
            
    else:
        f=plt.figure(1);
        f.subplots_adjust(hspace=0)
        f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=14)
        ax1 = plt.subplot2grid((2,2), (0, 0))
        ax2 = plt.subplot2grid((2,2), (1, 0))
        ax3 = plt.subplot2grid((2,2), (0, 1),rowspan=2)
        ax1.imshow(HDU_B.data[cy-rbkg/3.:cy+rbkg+1,cx-rbkg*.7:cx+rbkg*.7+1],origin='lower',vmin=np.min(bkg),vmax=np.max(bkg));
        ax2.imshow(-HDU_A.data[cy-rbkg:cy+rbkg/3+1.,cx-rbkg*.7:cx+rbkg*.7+1],origin='lower',vmin=np.min(bkg),vmax=np.max(bkg)); 
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
        P_nocor = AGPM.get_radial_profile(PSF_nocor, (xo,yo), nbin, disp=0)
        P_AGPM = AGPM.get_radial_profile(PSF_AGPM, (xo,yo), nbin, disp=0)
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

%pylab
# drawing sequence
#sequence = np.array([14,13,12,11,10,9,8,7,6,2,1,0]);
sequence = np.array([13,11,10,9,8,7,6,2,1,0]);
fsz = 9; #fonsize
#rim2 = int(2*FWHM);
vmin = 0.; vmax = 1.;
#vmin = -4;vmax = 0;  # for log10
lim = 15.*lbdDpp;

if obsSelect == 0:
    f=figure()
#    f.subplots_adjust(wspace=0)
    f.subplots_adjust(wspace=0.,right=.8);
    cbar_ax = f.add_axes([.81, 0.455, 0.007, 0.088]);
    #f.suptitle(target+', '+filt+' '+str(lb*1e6)+'um', fontsize=8)
    lamD = asecpmm/asecplbd*(114.50-np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05]))
    for k,i in enumerate(sequence): 
#        psf = PSFs[i,cim-rim2:cim+rim2+1,cim+offsets[i]-rim2:cim+offsets[i]+rim2+1]/np.max(PSFs);
        psf = ofaxMat[i]/np.max(ofaxMat)
#        psf = np.log10(psf); # for log10
        myAxes = [-lim+lamD[i],lim+lamD[i],-lim,lim];
        ax = subplot2grid((1,sequence.size),(0,k));
        im = ax.imshow(psf,origin='lower',vmin=vmin,vmax=vmax,extent=myAxes);
        ax.tick_params(labelsize=fsz,direction='out',length=3,right='off',top='off');
        ax.tick_params(left='off');
        ax.set_xticks([round(lamD[i],1)]);
        ax.set_yticks([-1,0,1]);  # is bugged, should appear on the figure...
#        ax.add_artist(Circle((lamD[i],0),0.5,color='w',fill=False));
        if k == 0:
            ax.tick_params(left='on');
            ax.set_ylabel(r"$\lambda$/D",fontsize=fsz)
            cb = f.colorbar(im, cax=cbar_ax);
            cb.ax.tick_params(labelsize=fsz);
            cb.set_ticks([vmin,(vmax+vmin)/2.,vmax]);
        if k == round(sequence.size/2):
            ax.set_xlabel(r"$\lambda$/D",fontsize=fsz);
#    text(-0.5,2.3,"log10",fontsize=fsz)
#    text(-15.,2.3,"circles = FWHM",fontsize=fsz)
    setp([a.get_yticklabels() for a in f.axes[1:]], visible=False);   
    savefig("CHLOUPY.png", dpi=300, transparent=True)


# UPDATE STOPS HERE !!!!
# ££££££££££££££££££££££

    
    f=plt.figure(2);
    xlabel(r"$\lambda$/D")
    ylabel("Radial profile")
    x1=asecpmm/asecplbd*(np.array([114.85,114.75,114.65,114.55,114.50,114.50,114.50,114.45,114.40,114.35,114.30,114.25,114.20,114.15,114.05])-114.50)
    y1=fluxMat/fluxNorm*0.838 # 0.838 = facteur d'Elsa
    x2 = np.array(sorted(abs(x1)))
    y2 = np.array(sorted(y1))
    dat = np.loadtxt('/Users/cdelacroix/Desktop/VISIR/VLT_circEP_Lyot_circ_aperture_nobsc_lyotout90_lyotin20_lbd37-1.txt')
    x3 = dat[:,0]
    y3 = dat[:,1]
    grid('on')
#    plot([-3,3], [.5,.5], 'k--')
    xlim([-3,3])
    plot(x3,y3, 'g', label="simulated")
    plot(x2,y2, '+m', markersize=8, mew=2, label="measured")
    yscale('log')
    plot([0,3], [.5,.5], 'k--')
    xlim([0,3])
    ylim(2e-2,1.5)
    
    IWA = 1.14 # 1.14 de Elsa
    plot([IWA,IWA], [2e-2,.5], 'k--', label='IWA = '+str(round(IWA,2))+'$\lambda$/D')
    xlabel(r"$\lambda$/D")
    legend(loc='lower right', fontsize=8)
    
    savefig(path+filt+'_offaxis.png', dpi=300, transparent=True);#plt.clf()


    
    




