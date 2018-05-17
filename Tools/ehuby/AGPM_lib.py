# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 11:37:06 2014

@author: ehuby
updates from cdelacroix
"""

import numpy as np
from astropy.io import fits         # fits file reading/writing
import matplotlib.pyplot as plt     
from matplotlib import cm           # palette for image display
import scipy.optimize as opt
import scipy.special as spe
import pylab


def get_image_VISIR(filename,i, obs, filt, nod, offaxis=0, disp=0):

    global target
#    obs=obsnum[i]
#    filt = filters[i]
#    nod = nodpos[i]
#    filename=glob.glob(dataDir+'*'+obs+'*')[0]
    """ Zone of interest """
    cx=390;cy=438; # obs 09-10,13-14
    if (obs=='0011') or (obs=='0012') or (obs=='0015') or (obs=='0016') or (obs=='0017') or (obs=='0018'):
        cy=550; # obs 11-12,15-18
    rbkg = 250
    rim = 100
    offset = -180*int((i%2==0)*2-1)
    """ Read image fits file """
    HDU=fits.open(filename)
    bkg = HDU[3].data[cy-rbkg:cy+rbkg,cx-rbkg:cx+rbkg]
    BCKm = np.median(bkg)
    PSF_AGPM = HDU[3].data[cy-rim:cy+rim,cx-rim:cx+rim] - BCKm
    PSF_nocor = -HDU[3].data[cy-rim+offset:cy+rim+offset,cx-rim:cx+rim] - BCKm
    HDU.close
    
    if disp == 1 :
        plt.imshow(bkg,origin='lower');plt.colorbar()
        plt.title('['+target+', '+filt+' filter, nod '+nod+'] Bkg median='+str(round(np.median(bkg),2)))
        plt.savefig(resDir+obs+'_bkg.png', dpi=300);plt.clf()
        plt.imshow(PSF_AGPM,origin='lower');plt.colorbar()
        plt.title('['+target+', '+filt+' filter, nod '+nod+'] AGPM centered')
        plt.savefig(resDir+obs+'_corON.png', dpi=300);plt.clf()
        plt.imshow(PSF_nocor,origin='lower');plt.colorbar()
        plt.title('['+target+', '+filt+' filter, nod '+nod+'] Off-axis PSF')
        plt.savefig(resDir+obs+'_corOFF.png', dpi=300);plt.clf()

    if offaxis==0:
        return PSF_AGPM
    elif offaxis==1:
        return PSF_nocor



def get_image_IRCT(filename, darkname, zoi=[0,-1,0,-1], disp=0, IRCT=1, rms=0): 
    ''' Computes dark subtracted median image contained in filename
        
        Read fits cube of images acquired on IRCT.

        [!] Data have to be read at the level 1 of HDU (not standard 0)
        
        Read fits cube of background images as well.
        
        Computes the median of PSF and background sequences.
        
        Returns the background subtracted median image.

        Optional key words:
            zoi:    zone of interest, defined by [x1,x2,y1,y2]
            
            disp:   display key word, 1 for displaying figures

            IRCT:   0 if not IRCT data (changes the level of HDU reading)

            rms:    display rms of the background images, for visual checking of bad pixels.    
    '''
    
    """ Zone of interest """
    x1=zoi[0]
    x2=zoi[1]
    y1=zoi[2]
    y2=zoi[3]
    
    """ Read image fits file """
    HDU=fits.open(filename)
    img=np.median(HDU[IRCT].data[:,y1:y2,x1:x2], axis=0)
    HDU.close

    """ Read dark fits file """
    HDU=fits.open(darkname)
    d=HDU[IRCT].data[:,y1:y2,x1:x2]
    N=(d.shape)[0]
    drk=np.median(d, axis=0)
    HDU.close
    
    """ Dark subtracted image """
    img_final=img-drk    
    
    if disp == 1 :
        plt.figure(num=1, figsize=(12,5))
        plt.clf()
        plt.subplot(231)
        plt.imshow(img)
        plt.colorbar()
        plt.title('Raw average image')
        plt.subplot(232)
        plt.imshow(drk)
        plt.colorbar()
        plt.title('Background image')
        plt.subplot(233)
        plt.imshow(img_final)
        plt.colorbar()
        plt.title('Final image')
    
    if rms == 1:
        # compute rms of the pixels to identify dead/hot picels
        rms=np.sqrt(np.mean(d**2, axis=0))
        rms_med=np.median(rms)
        plt.figure(2)
        plt.clf()
        plt.imshow(rms)#, vmin=0, vmax=rms_med*3)
        plt.colorbar()
        plt.title('RMS values over the '+str(N)+' images')
        print 'median rms = '+str(rms_med)
    
    return img_final
    

def get_image_YACA(allfilenames, zoi=[0,-1,0,-1], disp=0, rms=0): 
    '''Computes dark subtracted median image contained in filename
       
       Read fits cube of images acquired on YACADIRE.
       
       [!] Acquired data are already dark subtracted (no need for a "dark file")
       
       Computes the median of PSF sequence.
       
       Optional key words:
           zoi: zone of interest, defined by [x1,x2,y1,y2]
           
           disp:    display key word, 1 for displaying figures
           
           rms:     display rms of the PSF images (lack of background sequence), for visual checking of bad pixels.
    '''    
    
    """ Zone of interest """
    x1=zoi[0]
    x2=zoi[1]
    y1=zoi[2]
    y2=zoi[3]    
    
    """ Read image fits file """
    """ YACADIRE image: already dark subtracted """
    """ YACADIRE images are not cubes """
    N=len(allfilenames)
    
    if N==1 :
        HDU=fits.open(allfilenames[0])
        img_final=(HDU[0].data[x1:x2,y1:y2])
        HDU.close
    else :
        for k in range(0,N-1,1) :
            if k == 0 :
                HDU=fits.open(allfilenames[k])
                img=(HDU[0].data[x1:x2,y1:y2])[...,np.newaxis]
                HDU.close
            else :
                HDU=fits.open(allfilenames[k])
                d=HDU[0].data[x1:x2,y1:y2]
                img=np.concatenate((img,d[...,np.newaxis]), axis=2)
                HDU.close
        img_final=np.median(img,axis=2)
        #img_final=np.mean(img,axis=2)
    
    if disp == 1 :
        plt.figure(1)
        plt.clf()
        plt.imshow(img_final, cmap=cm.Greys_r, interpolation='none')
        plt.colorbar()
        plt.title('Final image')
        
    # compute rms of the pixels to identify dead/hot pixels
    if rms == 1:
        rms=np.sqrt(np.mean(img**2, axis=2))
        rms_med=np.median(rms)
        plt.figure(2)
        plt.clf()
        plt.imshow(rms, cmap=cm.Greys_r, interpolation='none', vmin=0, vmax=rms_med*3)
        plt.colorbar()
        plt.title('RMS values over the '+str(N)+' images')
        print 'median rms = '+str(rms_med)
    
    return img_final
    
def get_all_img_YACA(allfilenames, zoi=[0,-1,0,-1], disp=0, rms=0): 
    '''Computes dark subtracted median image contained in filename
       
       Read fits cube of images acquired on YACADIRE.
       
       [!] Acquired data are already dark subtracted (no need for a "dark file")
       
       Computes the median of PSF sequence.
       
       Optional key words:
           zoi: zone of interest, defined by [x1,x2,y1,y2]
           
           disp:    display key word, 1 for displaying figures
           
           rms:     display rms of the PSF images (lack of background sequence), for visual checking of bad pixels.
    '''    
    
    """ Zone of interest """
    x1=zoi[0]
    x2=zoi[1]
    y1=zoi[2]
    y2=zoi[3]    
    
    """ Read image fits file """
    """ YACADIRE image: already dark subtracted """
    """ YACADIRE images are not cubes """
    N=len(allfilenames)
    
    if N==1 :
        HDU=fits.open(allfilenames[0])
        img_final=(HDU[0].data[x1:x2,y1:y2])
        HDU.close
    else :
        for k in range(0,N-1,1) :
            if k == 0 :
                HDU=fits.open(allfilenames[k])
                img=(HDU[0].data[x1:x2,y1:y2])[...,np.newaxis]
                HDU.close
            else :
                HDU=fits.open(allfilenames[k])
                d=HDU[0].data[x1:x2,y1:y2]
                img=np.concatenate((img,d[...,np.newaxis]), axis=2)
                HDU.close
    
    return img
    

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    ''' Model function. 2D Gaussian.
    '''    
    
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
    
def oneD_Gaussian(x, amplitude, xo, sigma_x):
    ''' Model function. 1D Gaussian.
    '''    
    
    xo = float(xo)
    #g = offset + amplitude*np.exp( ((x-xo)/(2*sigma_x))**2 )
    g = amplitude*np.exp( -((x-xo)/(np.sqrt(2)*sigma_x))**2 )
    
    #print(amplitude, xo, sigma_x)
    
    return g
    
def poly6(x, x0, a0, a1, a2, a3, a4, a5, a6):
    ''' Model function. Polynomial function up to 6th order.
    '''
    
    xx = x-x0
    y = a0 + a1*xx + a2*xx**2 + a3*xx**3 + a4*xx**4 + a5*xx**5 + a6*xx**6
    
    return y

def poly6odd(x, x0, a0, a2, a4, a6):
    ''' Model function. Polynomial function up to 6th order.
    '''
    xx = x-x0
    y = a0 + a2*xx**2 + a4*xx**4 + a6*xx**6
    return y

def twoD_Airy((x,y), amplitude, xo, yo, F):
    ''' Model function. 2D Airy.
    '''    

    r = np.sqrt((x-xo)**2+(y-yo)**2)*F
    
    nx=r.shape[1]
    ny=r.shape[0]    
    maxmap=np.where(r==0, np.ones((ny,nx)), np.zeros((ny,nx)))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.unravel_index(maxmap.argmax(), maxmap.shape)
        r[indmax]=1.
    elif nbmax > 1:
        print 'ERROR in twoD_Airy: several nulls'
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1 :   
        Airy[indmax]=amplitude
    
    return Airy.ravel()

def oneD_Airy(x, amplitude, xo, F):
    ''' Model function. 1D Airy.
    '''

    r=(x-xo)*F
    nx=x.shape[0]
    
    maxmap=np.where(x==0, np.ones(nx), np.zeros(nx))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.argmax(maxmap)
        r[indmax]=1.
    elif nbmax > 1:
        print 'ERROR in oneD_Airy: several nulls'
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1 :   
        Airy[indmax]=amplitude
    
    return Airy
    
def oneD_Airy_log(x, amplitude, xo, F):
    ''' Model function. 1D log10(Airy).
    '''    

    r=(x-xo)*F
    nx=x.shape[0]
    
    maxmap=np.where(r==0, np.ones(nx), np.zeros(nx))
    nbmax=np.sum(maxmap)
    if nbmax == 1:
        indmax=np.argmax(maxmap)
        r[indmax]=1.
    elif nbmax > 1:
        print 'ERROR in oneD_Airy: several nulls'
    
    J=spe.jn(1, r)
    Airy=amplitude*(2*J/r)**2
    if nbmax == 1 :   
        Airy[indmax]=amplitude
    
    return np.log10(Airy)

def fit_gauss_2D(img):
    ''' Fits a 2D Gaussian pattern on the image.
        
        Returns the best fit parameters of the Gaussian shape.
        
        See twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''

    nx=img.shape[1]
    ny=img.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    x, y = np.meshgrid(x, y)

    init_xmax=np.unravel_index(img.argmax(), img.shape)[1]
    init_ymax=np.unravel_index(img.argmax(), img.shape)[0]
    initial_guess = (img.max(), init_xmax, init_ymax, 5, 5, 0, 0)
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), img.ravel(), 
                               p0=initial_guess)
    
    return popt

def fit_gauss_1D(y, x):
    ''' Fits a 1D Gaussian curve.
        
        Returns the best fit parameters of the Gaussian shape.
        
        See oneD_Gaussian(x, amplitude, xo, sigma_x, offset)
    '''

    #nx=y.shape[0]
    #x = np.linspace(0, nx-1, nx)

    init_xmax=x[y.argmax()]
    
    initial_guess = (y.max(), init_xmax, (x[-1]-x[0])/4.)
    popt, pcov = opt.curve_fit(oneD_Gaussian, x, y, p0=initial_guess)
    
    return popt
    
def fit_airy_2D(img, disp=0):
    ''' Fits a 2D Airy pattern on the image.
        
        Returns the best fit parameters of the Airy pattern.
        
        See twoD_Airy((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''

    nx=img.shape[1]
    ny=img.shape[0]
    x = np.linspace(0, nx-1, nx)
    y = np.linspace(0, ny-1, ny)
    x, y = np.meshgrid(x, y)

    init_xmax=np.unravel_index(img.argmax(), img.shape)[1]
    init_ymax=np.unravel_index(img.argmax(), img.shape)[0]
    initial_guess = (img.max(), init_xmax, init_ymax, .4)
    #initial_guess = (img[init_xmax, init_ymax]  , init_xmax, init_ymax, .6)

    #plt.figure(27)
    #plt.imshow(twoD_Airy((x,y), img[init_xmax, init_ymax]  , init_xmax, init_ymax, .6).reshape(nx,ny))

    popt, pcov = opt.curve_fit(twoD_Airy, (x, y), img.ravel(), p0=initial_guess)
    
    if disp != 0:
        data_fitted = twoD_Airy((x, y), *popt)
        #offset=np.min(img)
        plt.figure(disp)
        plt.clf()
        plt.subplot(1,3,1)
        plt.imshow(data_fitted.reshape(nx,ny), interpolation='none', cmap=cm.Greys_r)
        plt.colorbar()
        plt.subplot(1,3,2)
        plt.plot(img[popt[1],:])
        plt.plot(data_fitted.reshape(nx,ny)[popt[1],:], 'r--')
        plt.yscale('log')
        plt.subplot(1,3,3)
        plt.plot(img[:,popt[2]])
        plt.plot(data_fitted.reshape(nx,ny)[:,popt[2]], 'r--')
        plt.yscale('log')
        
        plt.figure(disp+1)
        plt.clf()
        plt.subplot(121)
        plt.imshow(data_fitted.reshape(nx,ny), interpolation='none', cmap=cm.Greys_r)
        plt.colorbar()      
        P_fit=get_radial_profile(data_fitted.reshape(nx,ny), (popt[1], popt[2]), 1, disp=10)
        P_mes=get_radial_profile(img, (popt[1], popt[2]), 1, disp=0)        
        plt.subplot(122)        
        plt.plot(P_mes, 'b')
        plt.plot(P_fit, 'r--')
        plt.yscale('log')
        
        print '\n--- Airy disk fit results ---'
        print 'Amplitude ='+str(popt[0])
        print 'Position of the maximum: \nxo='+str(popt[1])+' nyo='+str(popt[2])
        print 'F factor='+str(popt[3])
        print '-----------------------------'
    
    return popt
    
def fit_airy_1Dlog(Y, disp=0, initial_guess=[0.,0.,0.]):
    "fit one D   "
    
    nx=Y.shape[0]
    x=np.linspace(0,nx-1,nx)
    
#    minval=np.min(Y)
#    Y-=minval
#    ampl=np.max(Y)
#    Y=Y/ampl
    
    Ylog=np.log10(Y)
    
    if np.sum(initial_guess) == 0. :
        initial_guess=(np.max(Y), np.argmax(Y), .5)
    
    popt, pcov = opt.curve_fit(oneD_Airy_log, x, Ylog, p0=initial_guess, sigma=1./Y**2)  
    
    if disp != 0:
        data_fitted=oneD_Airy_log(x, *popt)        
        plt.figure(disp)
        plt.clf()
        plt.plot(Ylog, 'b')
        plt.plot(x, data_fitted, 'r--')
    
    return popt
    
def fit_airy_1D(Y, disp=0, initial_guess=[0.,0.,0.]):
    "fit one D   "
    
    nx=Y.shape[0]
    x=np.linspace(0,nx-1,nx)
    
#    minval=np.min(Y)
#    Y-=minval
#    ampl=np.max(Y)
#    Y=Y/ampl
    
    if np.sum(initial_guess) == 0. :
        initial_guess=(np.max(Y), np.argmax(Y), .5)
    
    popt, pcov = opt.curve_fit(oneD_Airy, x, Y, p0=initial_guess, sigma=1./Y)  
    
    if disp != 0:
        data_fitted=oneD_Airy(x, *popt)        
        plt.figure(disp)
        plt.clf()
        plt.plot(Y, 'b')
        plt.yscale('log')
        plt.plot(x, data_fitted, 'r--')
    
    return popt
    
def get_r_dist(nx,ny,xo,yo):
    ''' Returns the array of dimensions (nx,ny) with values corresponding the
        distance from the center (xo,yo). 
    '''
    
    x = np.linspace(0, nx, nx)-xo-1
    y = np.linspace(0, ny, ny)-yo-1
    x, y = np.meshgrid(x, y)
    
    return np.sqrt(x**2+y**2)
    
def get_radial_profile(img, (xo,yo), nbin, disp=0):
    ''' Computes the mean radial profile of the image.
    
        img:
            2D image.
        (xo,yo):
            center for the annuli.
        nbin:
            width of the annuli in pixels
        disp:
            optional key word for displaying the images.
            Its value will serve as the window number that will be created.
    '''
    
    (nx,ny)=img.shape
    r=get_r_dist(nx,ny,xo,yo)    
    
    r_max = np.max(r) # radius of the image
    r_max = np.max(r[xo,:])
    
    npts=int(r_max/nbin)
    O=np.ones((nx,ny))
    Z=np.zeros((nx,ny))
    Profile=np.zeros(npts)
    
    if disp != 0:
        plt.figure(disp)
        plt.clf()
        plt.subplot(121)
        plt.plot(xo,yo, 'xw')
        plt.title('PSF')
        plt.title('Averaged radial profile')
        plt.xlabel('Distance from center (pixels)')    
        plt.imshow(img, interpolation='none')
        plt.colorbar()
        
        val_min=np.min(img)
        val_max=np.max(img)        
        
        for k in range(0,npts-1,1):
            M=np.where(r>nbin*k, O, Z)*np.where(r<nbin*(k+1), O, Z)
            Profile[k]=np.sum(img*M)/np.sum(M)
            
            plt.figure(disp)
            plt.subplot(121)
            plt.imshow(img, interpolation='none')
            plt.pause(.005)
            plt.imshow(img*M, interpolation='none', vmin=val_min, vmax=val_max)
            plt.pause(.005)            
            plt.subplot(122)
            plt.plot(Profile, 'rx')
            plt.yscale('log')
        
        plt.plot(Profile, 'r')
        plt.subplot(121)
        plt.imshow(img, interpolation='none', vmin=val_min, vmax=val_max)
    
    for k in range(0,npts-1,1):
        M=np.where(r>nbin*k, O, Z)*np.where(r<nbin*(k+1), O, Z)
        Profile[k]=np.sum(img*M)/np.sum(M)   
    
    return Profile


def adjust_bckgr_level(img, xo, yo, R=0, disp=0):
    ''' Computes the median/mean background level of the image outside a 
        given radius.
        
        img:
            2D image
        (xo,yo):
            center of the PSF
        R: 
            radius of the circular zone to exclude.
        disp:
            optional keyword for displaying images.
    '''
    
    (nx,ny)=img.shape
    r=get_r_dist(nx,ny,xo,yo)
    
    M=np.double(img[np.where(r>R)])   
    
    bckgr_med=np.median(M)
    bckgr_mean=np.mean(M)
    
    
    print '\n----- Checking background level ------'
    print 'Background level = '+str(bckgr_med)+' [median]'
    print '                   '+str(bckgr_mean)+' [mean]'
    
    if disp != 0:
        plt.figure(disp)
        plt.clf()
        Pattern=np.where(r>R,np.ones((nx,ny)),np.zeros((nx,ny)))
        plt.imshow(Pattern*img, interpolation='none')
        plt.colorbar()
        plt.title('Background median on external area = '+str(bckgr_med))
    
    if bckgr_med != 0 :
        img=img-bckgr_med
        print ' -> background adjusted to 0 median'
    
    print '------------------------------------- '
        
    return (img, bckgr_med, bckgr_mean)
