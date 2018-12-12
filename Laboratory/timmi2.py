import matplotlib.pyplot as plt
import Tools.img_processing as impro
from astropy.io import fits
import os.path
import numpy as np

# some inputs
folder = os.path.expandvars('$HOME/INSTRUMENTS/TIMMI-2/')
names = ('near1', 'near2', 'agpmn3')
labels = {'near1':'AGPM-N-BT1', 'near2':'AGPM-N-BT2', 'agpmn3':'AGPM-N3'}
rim = 40
xbin = 1
ymin = 1e-6
pitch = 30*1e-6 # pixel pitch = 30 µm
#lam = 10.220*1e-6
#lam = 10.551*1e-6
lam = 10.675*1e-6
Fnum = 20

# some calculations
lamstr = str(int(lam*1e9))
pixscale = pitch/(lam*Fnum)
FWHM = 1/pixscale # Full Width at Half Maximum = lambda/D
(xo,yo) = (rim,rim)
(nx,ny) = (2*rim+1,2*rim+1)
circ_area = np.where(impro.get_r_dist(nx,ny,xo,yo) < FWHM/2) # radius of the integrated zone

# create null depth curves
x = np.arange(rim+1)*pixscale*xbin
for name in names:
    # file names
    file_noagpm = os.path.join(folder, name, '20181115_%smum_ls10mm_%sz925mm_no_01.fits'%(lamstr,name))
    file_offaxis = os.path.join(folder, name, '20181115_%smum_ls10mm_%sz925mm_offaxis_01.fits'%(lamstr,name))
    file_onaxis = os.path.join(folder, name, '20181115_%smum_ls10mm_%sz925mm_onaxis_01.fits'%(lamstr,name))

    # get data
    noagpm = -fits.getdata(file_noagpm)
    offaxis = -fits.getdata(file_offaxis)
    onaxis = -fits.getdata(file_onaxis)

    # subtract the background
    noagpm -= np.median(noagpm)
    offaxis -= np.median(offaxis)
    onaxis -= np.median(onaxis)

    # case 1: No AGPM
    res = impro.fit_gauss_2D(noagpm)
    (A1,yc,xc) = (res[0],int(res[1]),int(res[2]))
    noagpm_rim = noagpm[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
    y1 = impro.get_radial_profile(noagpm_rim, (xo,yo), xbin)
    I1 = np.sum(noagpm_rim[circ_area])

    # case 2: AGPM on-axis
    res = impro.fit_gauss_2D(offaxis)
    (A2,yc,xc) = (res[0],int(res[1]),int(res[2]))
    offaxis_rim = offaxis[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
    y2 = impro.get_radial_profile(offaxis_rim, (xo,yo), xbin)
    I2 = np.sum(offaxis_rim[circ_area])

    # case 3: AGPM off-axis
    res = impro.fit_gauss_2D(onaxis)
    (A3,yc,xc) = (res[0],int(res[1]),int(res[2]))
    onaxis_rim = onaxis[xc-rim+1:xc+rim,yc-rim+1:yc+rim]
    y3 = impro.get_radial_profile(onaxis_rim, (xo,yo), xbin)
    I3 = np.sum(onaxis_rim[circ_area])

    # normalize by the peak of the off-axis PSF
    y2 /= A2
    y3 /= A2

    # clip small/negative values to background level
    bkgd = np.var(onaxis/A2)
    y1[y1<bkgd] = np.min(y1[y1>ymin])
    y2[y2<bkgd] = np.min(y2[y2>ymin])
    y3[y3<bkgd] = np.min(y3[y3>ymin])

    # print stuff
    if True:
        print('1 = no agpm; 2 = off-axis psf; 3 = on-axis psf')
        print('Wavelength = %s nm'%lamstr)
        print('   %s Peak: (A1, A2, A3) = (%.1f, %.1f, %.1f)'%(name,A1,A2,A3))
        print('   %s Circ: (I1, I2, I3) = (%.1f, %.1f, %.1f)'%(name, I1,I2,I3))
        print('   %s Transmittance: A2/A1 = %.4f, I2/I1 = %.4f' %(name, A2/A1,I2/I1))
        print('   %s Rejection: A2/A3 = %d, I2/I3 = %d' %(name, A2/A3,I2/I3))

    # plot curve
    f2 = plt.figure(1)
    if name is 'near1':
        plt.plot(x, y2, 'k', label='off-axis PSF')
    else:
        plt.plot(x, y2, 'k')
    plt.plot(x, y3, label=labels[name])

"""figures"""
plt.yscale('log')
plt.xlim(0,6)
plt.ylim(ymin,1)
plt.plot([0.5,0.5],[ymin,1],'k:')
plt.xlabel('Angular separation in $x/(\lambda F\#)$')
plt.ylabel('Null Depth (normalized intensity)')
plt.title('Wavelength $\lambda = %s nm$'%lamstr)
plt.legend()
plt.grid()
plt.show(block=False)
plt.savefig(os.path.join(folder, '%snm.png'%lamstr), dpi=300, transparent=True)


#plt.figure(2)
#plt.imshow(noagpm, origin='lower')
#plt.colorbar()
#plt.show(block=False)
#plt.savefig(os.path.join(folder, 'img_noagpm.png'), dpi=300, transparent=True)

