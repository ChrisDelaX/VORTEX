import os, os.path
import numpy as np

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = 180 + np.rad2deg(np.arctan2(y, x))  # 0-360 [deg]
        return (rho, phi)

