from Tools.Misc import cart2pol
import numpy as np
import matplotlib.pyplot as plt

# Input paramters
# ---------------
N = 2**9 #10
#N2 = 2*N+1;
Rmax = 4 #2
D = 1 # lambda/D
lim = 8 #15
mask = 1 # 0=FQPM ; 1=AGPM
Lyot_rat = 0.9  # Lyot-stop ratio
Plim = 1.5
Flim = 3 #5#8#
Mlim = .8 #2
test = 1 # lambda/D
Rexp = 100 #


R = np.sqrt(N/np.pi) #32#1.5    # The outer radius -- the sqrt(N/pi) scale is roughly
                                # the best size for good sampling of both the aperture 
                                # and PSF
lamoverd = N/(2*R)              # number of samples in the image plane corresponding to
                                # lambda f-number

fontsize = 18
xy = np.arange(-N/2, N/2)
#xy = np.arange(-(N-1)/2., (N-1)/2.+1)
xyR = xy/R
xyF = xy/lamoverd
[X, Y] = np.meshgrid(xy, xy)
[RHO, TH] = cart2pol(X,Y)


""" Entrance pupil """

P = np.zeros((N,N))
P[np.where(RHO <= R)] = 1

plt.imshow(P,origin='lower',cmap="CMRmap");
plt.title('entrance pupil');
plt.xlabel('$D/2$');plt.ylabel('$D/2$')

P = exp(-(RHO/R).^Rexp);
LYOT = exp(-(RHO/(Lyot_rat*R)).^Rexp); % Lyot stop function

%modified pupil:
F = fftshift(fft2(fftshift(P)));
P2 = fftshift(fft2(fftshift(F.*M)));
P2(RHO <= R) = 0;
F2 = fftshift(ifft2(ifftshift(P2)));
PUPIL = fftshift(ifft2(ifftshift(F2./M)));
%clear P2 F2
%PUPIL = PUPIL./sqrt((max(max(abs(PUPIL).^2))));
    newFig;imagesc(xyR,xyR,abs(PUPIL));
    colormap(gray(256));colorbar;
    axis image;axis([-Plim Plim -Plim Plim]);
    title('{\bf entrance pupil}');xlabel('$D/2$');ylabel('$D/2$')
    %tick2latex;
    print('-depsc2',[folder '/' maskstr num2str(lp*mask) '_N=' num2str(N) '_R=' num2str(R,4) '_entrance_' folder '.eps'], '-r300');

