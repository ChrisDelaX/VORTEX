from Tools.Misc import cart2pol
import numpy as np
import matplotlib.pyplot as plt

""" Input parameters """

N = 2**10 #12
Cgoal = 1e-8                    # contrast goal
FPMradius = [3, 10.5]           # focal plane mask radius
FPMangle = 45                   # focal plane mask angle
Plim = 1.5                      # size of the pupil for display, in radii
Flim = 20                       # size of the image for display, in lam/D
Clim = -6                       # contrast limit for display
maxiter = 1000                   # maximum number of iterations

""" Sampling """

close('all')
R = 3*np.sqrt(N/np.pi)        # The outer radius -- the sqrt(N/pi) scale is roughly
                                # the best size for good sampling of both the aperture 
                                # and PSF.
lamD = N/(2*R)                  # number of samples in the image plane corresponding to
                                # lambda f-number

xy = np.arange(-N/2, N/2)
#xy = np.arange(-(N-1)/2., (N-1)/2.+1)
[X, Y] = np.meshgrid(xy, xy)
[RHO, TH] = cart2pol(X,Y)
AxesP = np.array([xy[0],xy[-1],xy[0],xy[-1]])/R
AxesF = np.array([xy[0],xy[-1],xy[0],xy[-1]])/lamD
RhoP = RHO/R
RhoF = RHO/lamD
xyP = xy/R
xyF = xy/lamD


""" Entrance pupil """

P = np.ones((N,N))
P[np.where(RhoP > 1)] = 0

plt.figure(1)
plt.imshow(P,origin='lower',interpolation='nearest',cmap="gray",extent=AxesP)
plt.colorbar()
plt.xlim(-Plim, Plim);plt.ylim(-Plim, Plim)
plt.title('entrance pupil')
plt.xlabel('$x/R$');plt.ylabel('$y/R$')


""" Focal plane """

F = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(P)))
PSF = np.log10(abs(F)/np.max(abs(F)))

plt.figure(2)
plt.imshow(PSF,origin='lower',interpolation='nearest',cmap="CMRmap",\
        extent=AxesF,vmin=Clim,vmax=0)
plt.colorbar()
plt.xlim(-Flim, Flim);plt.ylim(-Flim, Flim)
plt.title('point spread function')
plt.xlabel('$\lambda/D$');plt.ylabel('$\lambda/D$')


""" Contrast in the dark region """

mask = (RhoF > FPMradius[0]) * (RhoF < FPMradius[1]) * (abs(np.tan(TH)) < np.deg2rad(FPMangle))
C = np.mean(abs(F[mask])/np.max(abs(F)))


""" Mask optimization """

count = 0
plotfig = 0
while (C > Cgoal):# and count < maxiter:
    
    count += 1
    
    F[mask] = 0
#    F[~mask] = 0
    
    P = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(F)))
    P[np.where(RhoP > 1)] = 0
    P = abs(P)/np.max(abs(P))
    #P = np.round(P)
#    xk = sorted(range(len(P)**2), key=lambda k: P.reshape(len(P)**2)[k])
    
    
    if plotfig == 1:
        plt.figure(1)
        plt.imshow(P,origin='lower',interpolation='nearest',cmap="gray",extent=AxesP)
        plt.colorbar()
        plt.xlim(-Plim, Plim);plt.ylim(-Plim, Plim)
        plt.title('entrance pupil #'+str(count))
        plt.xlabel('$x/R$');plt.ylabel('$y/R$')
    
    F = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(P)))
    PSF = np.log10(abs(F)/np.max(abs(F)))
    
    if plotfig == 1:
        plt.figure(2)
        plt.imshow(PSF,origin='lower',interpolation='nearest',cmap="CMRmap",\
                extent=AxesF,vmin=Clim,vmax=0)
        plt.colorbar()
        plt.xlim(-Flim, Flim);plt.ylim(-Flim, Flim)
        plt.title('point spread function #'+str(count))
        plt.xlabel('$\lambda/D$');plt.ylabel('$\lambda/D$')
    
    C = np.mean(abs(F[mask])/np.max(abs(F)))
    print count, C


###

np.save('P.npy',P)
np.save('F.npy',F)

P = np.load('P.npy')
F = np.load('F.npy')

