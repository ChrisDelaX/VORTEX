import matplotlib.pyplot as plt
import Tools.img_processing as impro
from astropy.io import fits
import os.path
import numpy as np

""" some inputs """
folder = os.path.expandvars('$HOME/INSTRUMENTS/TIMMI2')
data = 'Data'               # name of the data folder
date = '20190201'           # measurement date
nsuffix = 2                 # number of suffix digits: 2 = '01', 4 = '0001',...
zfocus = 10                 # position of the AGPM: 10, 925,...
pitch = 30e-6               # pixel pitch in m
Fnum = 20                   # F number
xlim = (0, 6)               # figure x-axis limit: separation in x/(lam*f/D)
ylim = (1e-6, 1)            # figure y-axis limit: null depth
negative = True             # images are negative
print_stuff = True          # print values in the console
xbin = 1                    # radial profile binning parameter
rim = 25                    # cropped image radius in pixels
(xo,yo) = (rim,rim)
(nx,ny) = (2*rim+1,2*rim+1)

# select wavelengths
lams = ['10220', '10551', '10675'] 

# specifications
lams_specs = {'10220': {'names': ['BT2'],
                      'nimages': [5],
                   'last_deeps': [True]},
              '10551': {'names': ['BT1', 'BT2', 'BT3', 'N3'],
                      'nimages': [6, 6, 1, 6],
                   'last_deeps': [True, True, False, True]},
              '10675': {'names': ['BT2'],
                      'nimages': [6],
                   'last_deeps': [True]}}

labels = {'BT1': 'AGPM-N-BT1 (re-etched)',
          'BT2': 'AGPM-N-BT2 (re-etched)',
          'BT3': 'AGPM-N-BT3',
           'N3': 'AGPM-N3'}

colors = {'BT1': 'C0',
          'BT2': 'C1',
           'N3': 'C2',
          'BT3': 'C3'}
          
          
for lamstr in lams:
    
    # some calculations
    lam = float(lamstr)/1e9                                     # lambda in m
    pixscale = pitch/(lam*Fnum)                                 # pixel scale
    FWHM = 1/pixscale                                           # lambda/D
    circ_area = np.where(impro.get_r_dist(nx,ny,xo,yo) < FWHM/2)# integrated zone
    x = np.arange(rim+1)*pixscale*xbin                          # ang separations
    
    """ create null depth curves """
    lam_specs = lams_specs[lamstr]
    for iname, (name, nimage, last_deep) in enumerate(zip( \
            lam_specs['names'], lam_specs['nimages'], lam_specs['last_deeps'])):
        
        # initialize lists
        y1s, y2s, y3s = [], [], []
        Rs, Ts = [], []
        
        for i in range(nimage):
            
            # create name suffix
            suffix = ('%d'%(i+1)).zfill(nsuffix)
            # last image name ends with 'deep' (for a deep integration)
            if last_deep is True and i == nimage - 1:
                suffix += 'deep'
            
            # filenames
            nolamstr = lamstr#'10551'
            noagpm_name = '%s_%snm_ls10mm_%s_z%smm_noagpm_%s' \
                    %(date,nolamstr,name,zfocus,suffix)
            offaxis_name = '%s_%snm_ls10mm_%s_z%smm_offaxis_%s' \
                    %(date,lamstr,name,zfocus,suffix)
            onaxis_name = '%s_%snm_ls10mm_%s_z%smm_onaxis_%s' \
                    %(date,lamstr,name,zfocus,suffix)
            
            # load fits files from Data folder
            try:
                noagpm = fits.getdata(os.path.join(folder, data, '%s.fits'%noagpm_name))
            except:
                if negative is True:
                    noagpm = -noagpm
                pass # keep the last 'noagpm' image
            offaxis = fits.getdata(os.path.join(folder, data, '%s.fits'%offaxis_name))
            onaxis = fits.getdata(os.path.join(folder, data, '%s.fits'%onaxis_name))
            
            # change sign
            if negative is True:
                noagpm = -noagpm
                offaxis = -offaxis
                onaxis = -onaxis
            
            #find center of the off-axis image
            res = impro.fit_gauss_2D(offaxis)
            (xc,yc) = (int(res[2]),int(res[1]))
            
            # case 1: No AGPM
            noagpm_rim = noagpm[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
            noagpm_rim -= np.median(noagpm)
            y1 = impro.get_radial_profile(noagpm_rim, (xo,yo), xbin)
            A1 = np.max(y1)
            I1 = np.sum(noagpm_rim[circ_area])
            
            # case 2: AGPM off-axis
            offaxis_rim = offaxis[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
            offaxis_rim -= np.median(offaxis)
            y2 = impro.get_radial_profile(offaxis_rim, (xo,yo), xbin)
            A2 = np.max(y2)
            I2 = np.sum(offaxis_rim[circ_area])
            
            # case 3: AGPM on-axis
            onaxis_rim = onaxis[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
            onaxis_rim -= np.median(onaxis)
            y3 = impro.get_radial_profile(onaxis_rim, (xo,yo), xbin)
            A3 = np.max(y3)
            I3 = np.sum(onaxis_rim[circ_area])
            
            # plot PSFs
            if True:
                fig, (ax1, ax2) = plt.subplots(1,2)
                im1 = ax1.imshow(np.clip(np.log(offaxis_rim/A2), -7, 0), \
                        origin='lower', cmap="CMRmap")
                im2 = ax2.imshow(np.clip(np.log(onaxis_rim/A2), -7, 0), \
                        origin='lower', cmap="CMRmap")
                ax1.set_title('FWHM=%d'%(I2))
                ax2.set_title('R=%d'%(I2/I3))
                fig.subplots_adjust(right=0.8)
                cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                cb = fig.colorbar(im1, cax=cbar_ax)
                plt.show(block=False)
                plt.savefig(os.path.join(folder, '%s.png'%offaxis_name), \
                        dpi=300, transparent=True)
                plt.close()
            
            # normalize by the peak of the off-axis PSF
            y2 /= A2
            y3 /= A2
            
            # calculate rejection and transmittance
            Rs.append(I2/I3)
            Ts.append(I2/I1)
            
            # clip small/negative values to background level
            bkgd = np.var(onaxis/A2)
            y1[y1<ylim[0]] = np.min(y1[y1>bkgd])
            y2[y2<ylim[0]] = np.min(y2[y2>bkgd])
            y3[y3<ylim[0]] = np.min(y3[y3>bkgd])
            
            # append curve to list
            y1s.append(y1)
            y2s.append(y2)
            y3s.append(y3)
            
            # print stuff
            if print_stuff is True:
                print('1 = no agpm; 2 = off-axis psf; 3 = on-axis psf')
                print('Wavelength = %s nm'%lamstr)
                print('   %s Peak: (A1, A2, A3) = (%.1f, %.1f, %.1f)'%(name,A1,A2,A3))
                print('   %s Circ: (I1, I2, I3) = (%.1f, %.1f, %.1f)'%(name, I1,I2,I3))
                print('   %s Transmittance: A2/A1 = %.4f, I2/I1 = %.4f' %(name, A2/A1,I2/I1))
                print('   %s Rejection: A2/A3 = %d, I2/I3 = %d' %(name, A2/A3,I2/I3))
            
            # plot all curves in dotted lines
            plt.figure(1)
            plt.plot(x, y2, 'k:')
            plt.plot(x, y3, ':', color=colors[name])
        
        # plot mean curves
        if iname == 0:
            plt.plot(x, np.median(y2s,0), 'k', label='off-axis PSF')
        else:
            plt.plot(x, np.median(y2s,0), 'k')
        plt.plot(x, np.median(y3s,0), color=colors[name], label='%s: R=%d, T=%.2f' \
    #            %(labels[name], np.median(Rs), np.nan))
                %(labels[name], np.median(Rs), np.median(Ts)))
    
    """figures"""
    plt.yscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.plot([0.5,0.5],ylim,'k:') # FWHM
    plt.xlabel('Angular separation in $x/(\lambda F\#)$')
    plt.ylabel('Null Depth (normalized intensity)')
    plt.title('Wavelength $\lambda = %s nm$'%lamstr)
    plt.legend()
    plt.grid()
    plt.show(block=False)
    plt.savefig(os.path.join(folder, '%s_%snm_.png'%(date, lamstr)), dpi=300, transparent=True)
    plt.close()

