# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:48:45 2015

Comoving Box Visualization
@author: Sheridan Beckwith Green (sbg@unc.edu)
"""

import math
import numpy as np
import sys

def distance(x,y,z):
    return np.sqrt(x*x + y*y + z*z)

#parameters
#make this parameter one that we load in from command line or parameter file
vis_size = 256 #the dimension of our cube for visualization (could eventually make a parameters file, this would go in)

#read in our file
infile    = open(sys.argv[1],"rb")
latticefile   = open(sys.argv[1]+".vis","w")
particlefile = open(sys.argv[1]+".particle","w")

[xcenter,ycenter,zcenter] = sys.argv[2:5] #halo center
rad = float(sys.argv[5]) #halo radius

#### Begin header fortran read statement
dummy = np.fromfile(infile,dtype=np.int32,count=1)

# The following should add up to 256 bytes
npart     = np.fromfile(infile,dtype=np.int32,count=6)
mass      = np.fromfile(infile,dtype=np.float64,count=6)
a         = np.fromfile(infile,dtype=np.float64,count=1)
z         = np.fromfile(infile,dtype=np.float64,count=1)
flagsfr   = np.fromfile(infile,dtype=np.int32,count=1)
flagfeedb = np.fromfile(infile,dtype=np.int32,count=1)
nall      = np.fromfile(infile,dtype=np.int32,count=6)
useless   = np.fromfile(infile,dtype=np.int32,count=2)
sim_size  = np.fromfile(infile,dtype=np.float64,count=1)*1000 #in pc/h
omega0    = np.fromfile(infile,dtype=np.float64,count=1)
omega_l   = np.fromfile(infile,dtype=np.float64,count=1)
hubble    = np.fromfile(infile,dtype=np.float64,count=1)
unused    = np.fromfile(infile,dtype=np.int32,count=24)

#### End header fortran read statement
dummy = np.fromfile(infile,dtype=np.int32,count=1)

#### Begin positions fortran read statement
dummy = np.fromfile(infile,dtype=np.int32,count=1)

# Here we will assume that there are only particles of type 1 (x1000 to convert to parsecs/h)
positions = np.fromfile(infile,dtype=np.float32,count=npart[1]*3).reshape([npart[1],3])*1000
#Positions of the particles are in COMOVING, as well as in pc/h

#move positions to be relative to our radius
#these are the raw positions
positions[:,0] = positions[:,0] - float(xcenter)
positions[:,1] = positions[:,1] - float(ycenter)
positions[:,2] = positions[:,2] - float(zcenter)
#may need to rethink this when we're generating a density profile instead of just visualization.

#now we have positions, we need to get the box dimensions and split it into a course mesh
res = rad / vis_size #this is the size of one "pixel"

#positions binned into the "pixels"
pos_bins = positions / res #when used as index, numpy takes the floor of these values 
cube = np.ndarray((vis_size,vis_size,vis_size)) #goes from 0 to 255

print "Generating 3d lattice, placing particles in lattice..."
for pos in positions:
    #if(distance(pos[0],pos[1],pos[2]) < rad):
    if((np.abs(pos[0]) < rad) & (np.abs(pos[1]) < rad) & (np.abs(pos[2]) < rad)):
        xbin = ((pos[0] / res) + vis_size) / 2.0 #defines the center of our box to be at 0 radius
        ybin = ((pos[1] / res) + vis_size) / 2.0
        zbin = ((pos[2] / res) + vis_size) / 2.0
        cube[xbin][ybin][zbin] += 1
        particlefile.write(str(pos[0])+" " +str(pos[1])+" " +str(pos[2])+"\n")
    
#need to store the location of particle as well for data analysis purposes

print "Writing to file"
#we can normalize at end in visualization script
for i in range(0,vis_size):
    for j in range(0,vis_size):
        for k in range(0,vis_size):
            latticefile.write(str(cube[i][j][k])+"\n")

