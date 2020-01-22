import numpy as np
from numpy.linalg import norm


# map of size 2*N + 1
N = int(100)
xy = np.array(range(-N,N+1))

# 2D map: cartesian and polar coordinates
X, Y = np.meshgrid(xy,xy)
THETA = np.arctan2(Y, X)
RHO = np.sqrt(X**2 + Y**2)

# optical axis
OA = (THETA + np.pi/2) % np.pi - np.pi/2
# calculate TE and TM, for input polarization ux = 1, uy = 0
TM = np.array([np.cos(OA),np.sin(OA)]) * np.cos(OA)
TE = np.array([np.sin(OA),-np.cos(OA)]) * np.sin(OA)

# intensity and phase mismatch
tetm_rat = 1.2 # TM/TE
tetm_phi = 0.14*np.pi
# rotation matrix
rot = np.array([[np.cos(tetm_phi), -np.sin(tetm_phi)], \
        [np.sin(tetm_phi), np.cos(tetm_phi)]])

# reconstruct the output polarization
angles = np.zeros((2*N+1,2*N+1))
for i in xrange(2*N+1):
    for j in xrange(2*N+1):
        output = TM[:,i,j] + np.dot(rot,TE[:,i,j]/tetm_rat)
        angles[i,j] = np.arctan2(output[1],output[0])

# figure
%pylab
myAxes = [-N,N,-N,N]
imshow(angles,origin='lower',extent=myAxes,cmap='hsv')#,clim=[-pi,pi])
colorbar()
