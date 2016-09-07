# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:39:55 2015

@author: greensb
"""

import numpy as np
import sys
from mayavi import mlab

vis_size = 256 #the dimension of our cube for visualization
cube = np.ndarray((vis_size,vis_size)) #goes from 0 to 255, initialized to zeros

#generate flat image from final data
#figure out way to pick which dimension to look at
#ALSO, NEED A WAY TO DETERMINE THE VISUALIZATION SIZE BASED ON PA

maxval = 0
data = open(sys.argv[1],"rb")
for j in range(0,vis_size):
    for k in range(0,vis_size):
        cube[j][k] = 0 #initialize it to zero
        for m in range(0,vis_size):
            cube[j][k] = cube[j][k] + float(data.readline())
        if(cube[j][k] > maxval):
            maxval = cube[j][k]

print maxval
cube = (cube / maxval) * 255
            
#mlab.pipeline.volume(mlab.pipeline.scalar_field(cube),vmin=0.2, vmax=0.8)
mlab.imshow(cube,vmin=0.0,vmax=180.0)
mlab.show()