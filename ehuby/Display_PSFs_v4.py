# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 09:04:26 2014

@author: ehuby
updates from cdelacroix

v2: can be used for YACADIRE & IRCT data
v3: backgroun subtraction for 2015-01-09 YACADIRE data in K-band
v4: for VISIR observations (added by cdelacroix)
"""

import numpy as np
import glob
from astropy.io import fits
#from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.optimize as opt
import pylab as pl
import scipy.special as spe
import time

#from scipy.fftpack import fft2, ifft2

import AGPM_lib as AGPM

num_back = []
obsnum =['']
num_nocor =['']

from all_PSFs_VISIR import *



save_png=1

#prefix=date+'_'+name_AGPM+'_'+face+'_'+source+'_allPSF'
today=time.strftime("%Y%m%d")
#dirSAVE='/home/ehuby/VORTEX/Results/Tests_AGPM/'+name_AGPM+'/null_estimation_'+source+'/'+prefix

# constrained minimal value for the radial values, to allow log scale
minprofiles = 1e-6 
# binning parameter for the radial profile
nbin=1.



""" ******************** LOADING non attenuated PSFs ********************** """
""" *********************************************************************** """
k=0
file_name = glob.glob(dataDir+'*'+obsnum[k]+'*')[0]
PSF_AGPM = AGPM.get_image_VISIR(file_name,k,obsnum[k],filters[k],nodpos[k],0)
PSF_nocor = AGPM.get_image_VISIR(file_name,k,obsnum[k],filters[k],nodpos[k],1)
(nx,ny) = PSF_nocor.shape  
  


""" **** Fitting an Airy disk to find the coordinates of the center ******* """
# the fit is performed on the non attenuated PSF
# first guess for the coordinates of the center
xinit=np.unravel_index(PSF_nocor.argmax(), PSF_nocor.shape)[1]
yinit=np.unravel_index(PSF_nocor.argmax(), PSF_nocor.shape)[0]
# fitting routine
(A,xo,yo,Fairy)=AGPM.fit_airy_2D(PSF_nocor, disp=0)

# Note: the fit is not really efficient to estimate the width of the PSF
# I tunned this parameter manually
Fairy=3.4               # HeNe on IRCT  
Rairy=3.83/Fairy        # First zero of the Airy pattern in pixels
Rairy2=7.00/Fairy       # Second zero
Rairy3=10.17/Fairy      # Third zero
FWHM=Rairy/1.22         # Full Width at Half Maximum = lambda/D
lbdDpp=1/(Rairy/1.22)   # conversion factor in lbd/D per pixel

""" ****** NULL DEPTH ESTIMATION FROM PHOTOMETRY IN A CIRCULAR ZONE ******* """ 
""" *********************************************************************** """
Rcirc=FWHM/2.           # radius of the integrated zone
r=AGPM.get_r_dist(nx,ny,xo,yo)
circ_area = np.where(r<Rcirc)
Flux_nocor_circ=np.sum(PSF_nocor[circ_area])
Flux_AGPM_circ=np.sum(PSF_AGPM[circ_area])
N=round(Flux_nocor_circ/Flux_AGPM_circ,1)
print('Null depth (aperture photometry) = '+str(N))
""" *********************************************************************** """


# PSF display
plt.figure(10)
plt.clf()
plt.imshow(np.log10(PSF_nocor-np.min(PSF_nocor)),origin='lower')
plt.title("Mean off-axis PSF")
cb=plt.colorbar()
cb.set_label('log10')

# PSF display
plt.figure(20)
plt.clf()
plt.imshow(np.log10(PSF_AGPM-np.min(PSF_AGPM)),origin='lower')
plt.title("Mean attenuated PSF")
cb=plt.colorbar()
cb.set_label('log10')


""" *************************** RADIAL PROFILES *************************** """
""" *********************************************************************** """
# Note that these profiles are for display purposes only.
# the null depth is estimated from the photometry directly in the image
# there are two estimations of the mean radial profile (consistency check)
#    - 1) radial profile of the mean PSF
#    - 2) mean of all the individual profiles


""" 1) RADIAL PROFILE OF THE MEAN PSF """
# radial profile from mean PSF
P_nocor=AGPM.get_radial_profile(PSF_nocor, (xo,yo), nbin, disp=0)
P_AGPM=AGPM.get_radial_profile(PSF_AGPM, (xo,yo), nbin, disp=0)

# normalization:
norm=np.max(P_nocor)
P_nocor = P_nocor/norm
P_AGPM  = P_AGPM/norm

# check if the minimal value is negative (to allow correct log display)
# if this is the case, the minimal value is shifted to 'minprofile' value
valmin=np.min(P_AGPM)-minprofiles
if valmin < 0 :
    P_AGPM = P_AGPM - valmin

# corresponding radial distance
(nxp)=P_nocor.shape[0]
Xlbd=np.linspace(0,nxp-1,nxp)*lbdDpp

if bench!='VISIR':
    """ 2) MEAN RADIAL PROFILE from all PSF """
    """ non attenuated profiles """
    Allpro_nocor=np.zeros((nxp, N_nocor))
    
    plt.figure(16)
    plt.clf()
    for i in range(N_nocor):
        Allpro_nocor[:,i]=AGPM.get_radial_profile(AllPSF_nocor[:,:,i]-BCKm, (xo,yo), nbin, disp=0)
        Allpro_nocor[:,i]=Allpro_nocor[:,i]/norm
        
        # check if the minimal value is negative (to allow correct log display)
        # if this is the case, the minimal value is shifted to 'minprofile' value
        valmin = np.min(Allpro_nocor[:,i]) - minprofiles
        if valmin < 0:
            Allpro_nocor[:,i]=Allpro_nocor[:,i] - valmin + minprofiles

        # display
        plt.figure(16)
        plt.plot(Allpro_nocor[:,i], 'r--')
    
    
    """ attenuated profiles """
    Allpro_AGPM=np.zeros((nxp, N_AGPM))
    
    plt.figure(15)
    plt.clf()
    for i in range(N_AGPM):
        Allpro_AGPM[:,i]=AGPM.get_radial_profile(AllPSF_AGPM[:,:,i]-BCKm, (xo,yo), nbin, disp=0)
        Allpro_AGPM[:,i] = Allpro_AGPM[:,i]/norm
        
        # check if the minimal value is negative (to allow correct log display)
        # if this is the case, the minimal value is shifted to 'minprofile' value
        valmin = np.min(Allpro_AGPM[:,i])
        if valmin < 0 :
            Allpro_AGPM[:,i]=Allpro_AGPM[:,i] - valmin + minprofiles
        
        # display
        plt.figure(15)
        plt.plot(Allpro_AGPM[:,i], 'b--')
    
    
    P_nocor_mean = np.mean(Allpro_nocor, 1)
    P_AGPM_mean = np.mean(Allpro_AGPM, 1)
else:
    P_nocor_mean = P_nocor
    P_AGPM_mean = P_AGPM

""" display non attenuated profiles """
plt.figure(16)
plt.plot(P_nocor, 'r')         # radial profile of the mean PSF
plt.plot(P_nocor_mean, ':r')   # mean radial profile from all profiles
plt.title(name_AGPM+' '+face+' '+source+" - "+date+" - All off-axis radial profiles")
plt.yscale('log')
plt.ylim(1e-6, 1.)
""" display attenuated profiles """
plt.figure(15)
plt.plot(P_AGPM, 'b')         # radial profile of the mean PSF
plt.plot(P_AGPM_mean, ':b')   # mean radial profile from all profiles
plt.title(name_AGPM+' '+face+' '+source+" - "+date+" - All attenuated radial profiles")
plt.yscale('log')
plt.ylim(-1e-6,1)
""" display both mean profiles """
plt.figure(12)
plt.clf()
plt.plot(Xlbd*nbin, P_nocor, 'r', label="Off-axis PSF")
plt.plot(Xlbd*nbin, P_AGPM, 'b', label="Attenuated PSF (R="+str(int(N))+")")
plt.xlabel(r"$\lambda$/D")
plt.ylabel("Normalized PSF radial profile")
plt.legend(loc='upper right')
plt.title(name_AGPM+' '+face+' '+source+" - "+date)
plt.yscale('log')
plt.ylim(1e-6,1.)
plt.xlim(0,6.)
# markers for the PSF zeros
R=Rcirc*lbdDpp
plt.plot([R,R], [10**(-6), 1.], 'k--')
R=Rairy*lbdDpp
plt.plot([R,R], [10**(-6), 1.], 'm:')
R=Rairy2*lbdDpp
plt.plot([R,R], [10**(-6), 1.], 'm:')
R=Rairy3*lbdDpp
plt.plot([R,R], [10**(-6), 1.], 'm:')


if save_png == 1:
    plt.figure(12)
    plt.savefig(resDir+'_Radial_profiles.png', dpi=300)
    plt.figure(20)
    plt.savefig(resDir+'_mean_PSF_AGPM.png', dpi=300)
    plt.figure(10)
    plt.savefig(resDir+'_mean_PSF_nocor.png', dpi=300)    
    plt.figure(15)
    plt.savefig(resDir+'_all_profiles_AGPM.png', dpi=300)   
    plt.figure(16)
    plt.savefig(resDir+'_all_profiles_nocor.png', dpi=300)
