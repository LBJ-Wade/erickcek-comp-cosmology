# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 14:57:52 2015

@author: greensb
"""

import numpy as np
import sys
from mayavi import mlab

vis_size = 257 #the dimension of our cube for visualization
#cube = np.ndarray((vis_size,vis_size,vis_size)) #goes from 0 to 255, initialized to zeros
cube = np.sum(np.load(sys.argv[1]),axis=0)

#maxval = 0
maxval = np.max(cube)
#data = open(sys.argv[1],"rb")
#for j in range(0,vis_size):
#    for k in range(0,vis_size):
#        for m in range(0,vis_size):
#            cube[j][k][m] = float(data.readline())
#           if(cube[j][k][m] > maxval):
#                maxval = cube[j][k][m]

cube = (cube / maxval) * 255
            
mlab.pipeline.volume(mlab.pipeline.scalar_field(cube),vmin=0.2, vmax=0.8)
mlab.show()