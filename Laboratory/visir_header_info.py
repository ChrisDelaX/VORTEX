from Tools.ObservingBlock import ObservingBlock

folder = '/Users/cdelacroix/INSTRUMENTS/VISIR/2015-03-09_off-axis'
OB = ObservingBlock(folder, start='corono')
OB.__dict__.keys()
#print OB

# ERIS pupil
folder = '/Users/cdelacroix/INSTRUMENTS/VISIR'
OB = ObservingBlock(folder, start='pup')
%pylab
imshow(OB.getData(0))


#OB.getAttribute('ORIGFILE') # 'corono_0019.fits', 'corono_0020.fits', 'corono_0021.fits',
#        'corono_0022.fits', 'corono_0023.fits', 'corono_0024.fits',
#        'corono_0025.fits', 'corono_0026.fits', 'corono_0027.fits',
#        'corono_0028.fits', 'corono_0029.fits', 'corono_0030.fits',
#        'corono_0031.fits', 'corono_0032.fits', 'corono_0033.fits',
#        'corono_0034.fits', 'corono_0035.fits', 'corono_0036.fits',
#        'corono_0037.fits', 'corono_0038.fits', 'corono_0039.fits',
#        'corono_0040.fits', 'corono_0041.fits', 'corono_0042.fits',
#        'corono_0043.fits', 'corono_0044.fits', 'corono_0045.fits',
#        'corono_0046.fits', 'corono_0047.fits', 'corono_0048.fits',
#        'corono_0049.fits', 'corono_0050.fits', 'corono_0051.fits',
#        'corono_0052.fits', 'corono_0053.fits', 'corono_0054.fits'
#OB.getAttribute('EXPTIME') # 0.0024887
#OB.getAttribute('HIERARCH ESO DET CHOP TIM') # 0.0024887
#OB.getAttribute('HIERARCH ESO DET CHOP TRANSTIM') # 0.025
#OB.getAttribute('HIERARCH ESO DET SEQ1 WIN STRX') # 1 / Start-X (lower left)
#OB.getAttribute('HIERARCH ESO DET SEQ1 WIN STRY') # 388 / Start-Y (lower left)
#OB.getAttribute('HIERARCH ESO DET SEQ1 WIN NX') # 1024 / number of pixels along x
#OB.getAttribute('HIERARCH ESO DET SEQ1 WIN NY') # 250 / number of pixels along y
#OB.getAttribute('HIERARCH ESO DET SEQ1 REALDIT') # 0.0024887 / Integration time
#OB.getAttribute('HIERARCH ESO DET SEQ1 EXPTIME') # [12.4460887, 12.4460887, 12.4460887, 
#       12.4460887, 12.4460887, 12.4460887, 12.4460887, 12.4460887, 12.4460887, 12.4460887, 
#       24.8896887, 24.8896887, 24.8896887, 37.3332887, 37.3332887, 74.6640887, 74.6640887, 
#       124.4384887, 124.4384887, 124.4384887, 124.4384887, 49.7768887, 49.7768887, 24.8896887, 
#       24.8896887, 24.8896887, 24.8896887, 24.8896887, 24.8896887, 24.8896887, 24.8896887, 
#       24.8896887, 24.8896887, 49.7768887, 49.7768887, 49.7768887]

#OB.getAttribute('HIERARCH ESO DET NDIT') # 5000,  5000,  5000,  5000,  5000,  5000,  5000,  5000,  5000,
#        5000, 10000, 10000, 10000, 15000, 15000, 30000, 30000, 50000,
#       50000, 50000, 50000, 20000, 20000, 10000, 10000, 10000, 10000,
#       10000, 10000, 10000, 10000, 10000, 10000, 20000, 20000, 20000

#OB.getAttribute("CRPIX1", 1) # 0.5 / Ref pixel in axis1
#OB.getAttribute("CRPIX2", 1) # -386.5 / Ref pixel in axis2
#OB.getAttribute("HIERARCH ESO DET CHIP NAME", 1) # 'Aquarius' / Detector chip name
#OB.getAttribute("HIERARCH ESO DET CHIP DATE", 1) # '2012-06-20' / Date of installation
#OB.getAttribute("HIERARCH ESO DET CHIP PXSPACE", 1) # 3.000E-05 / Pixel-Pixel Spacing
#OB.getAttribute("HIERARCH ESO DET CHIP GAIN", 1) # 20.00 / Gain in e-/ADU
#OB.getAttribute("HIERARCH ESO DET CHIP RON", 1) # 300 / Read-out noise in e-




