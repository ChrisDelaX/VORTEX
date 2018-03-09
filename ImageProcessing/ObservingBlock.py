import os, os.path
import numpy as np

import astropy.units as u
from astropy.io import fits




class ObservingBlock(object):

    def __init__(self, folder, start=''):
        
        self.name = os.path.split(folder)[-1]
        self.start = start
        self.files = np.array([f for f in os.listdir(folder) if f[-5:] == '.fits' \
                and f[:len(start)] == start])
        
#        with fits.open(filename) as hdulist:
        
        return

    def __str__(self):
        """Use the command 'print' to print the attribute values contained in the object"""
        
        for att in self.__dict__.keys():
            print('%s: %r' % (att, getattr(self, att)))
        
        return 'ObservingBlock class object attributes'