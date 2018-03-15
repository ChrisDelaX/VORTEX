import os, os.path
import numpy as np
import glob as glob
from astropy.io import fits




class ObservingBlock(object):

    def __init__(self, folder, start='', seq=None):
        
        # file names
        if seq is not None:
            seq = np.array(seq, ndmin=1)
            self.files = np.array([glob.glob(os.path.join(folder, '*'+seq[i]+'*.fits'))[0] \
                    for i in range(seq.size)])
        else:
            self.files = np.array([f for f in os.listdir(folder) if f[-5:] == '.fits' \
                    and f[:len(start)] == start])
        
        # copy attributes
        self.folder = folder
        self.start = start
        self.seq = seq
        self.nfiles = self.files.size
            
        # corresponding HDU lists (does not store data, only headers)
        self.hdulists = np.empty(self.nfiles, dtype=object)
        for i,filename in enumerate(self.files):
            with fits.open(os.path.join(folder, filename)) as hdulist:
                self.hdulists[i] = hdulist
        
        return

    def __str__(self):
        """Use the command 'print' to print the attribute values contained in the object"""
        
        for att in self.__dict__.keys():
            print('%s: %r' % (att, getattr(self, att)))
        
        return 'ObservingBlock class object attributes'

    def getAttribute(self, att, ind=0):
        
        res = np.array([h[ind].header.get(att) for h in self.hdulists])
        
        return res
        
