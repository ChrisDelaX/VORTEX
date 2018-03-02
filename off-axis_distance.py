# -*- coding: utf-8 -*-
"""
Created on Mon May 11 07:20:50 2015

@author: cdelacroix
"""


import numpy as np
import AGPM_lib_visir as AGPM
import matplotlib.pyplot as plt  

dist = [114.05,114.15,114.2,114.25,114.3,114.35,114.4,114.65,114.75,114.85]

x = [373.39,377.54,379.64,381.57,383.32,384.72,384.46,400.18,403.58,407.72]
y = [164.43,164.45,164.53,164.55,164.61,164.58,164.67,164.3,163.72,163.51]
ymed = np.median(y[:7])
xOffset = ['114.85','114.75','114.65','114.55','114.50','114.50','114.50','114.45','114.40','114.35','114.30','114.25','114.20','114.15','114.05']    
dist2 = np.zeros(np.size(xOffset))
x2 = np.zeros(np.size(xOffset))
y2 = np.zeros(np.size(xOffset))
for k in range(np.size(xOffset)):
    dist2[k] = float(xOffset[k])
    x2[k] = 372.89+.429125*(1+(int(xOffset[k][-2:])-5))
    y2[k] = 164.55 #ymed    

plt.figure(1); plt.clf()
plt.plot([dist[0],dist[-1]],[x[0],x[-1]])
plt.plot(dist[:7],x[:7], 'r+-', label="Off-axis X-distance")
plt.plot(dist[7:],x[7:], 'r+-', label="Off-axis X-distance")
plt.plot(dist2,x2,'g+')
plt.xlabel("x offset (mm)"); plt.ylabel("x position (pixels)")
plt.savefig('xo-position.png', dpi=300);plt.clf()


plt.figure(2); plt.clf()
#plt.plot([dist[0],dist[-1]],[y[0],y[0]])
plt.plot([dist[0],dist[-1]],[y2[0],y2[0]])
plt.plot(dist[:7],y[:7], 'r+-', label="Off-axis X-distance")
plt.plot(dist[7:],y[7:], 'r+-', label="Off-axis X-distance")
plt.ylim(160,170)
plt.plot(dist2,y2,'g+')
plt.xlabel("x offset (mm)"); plt.ylabel("y position (pixels)")
plt.savefig('yo-position.png', dpi=300);plt.clf()

