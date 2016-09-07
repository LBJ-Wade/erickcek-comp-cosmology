# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:58:33 2015

@author: greensb
PROCEDURE FOR VISUALIZING MICROHALOS
1. Run this script, giving it as arguments the corresponding ROCKSTAR halo file and the associated gadget snapshot blocks
2. Using the *.vis files, run block_combine.py to combine them all into one .txt file for visualization.
3. Use flat_visualize.py and visualize.py with the .txt file to view microhalo

PROCEDURE FOR GENERATING DENSITY PROFILES
1. Run this script, same as above, and it will output the locations of the particles (we still need to get particle mass somehow)
2. Run command cat *.particle > all_particles.txt
3. Run density_prof.py, which takes all the particles, bins them, and outputs a density profile file, which we can plot on a personal computer
"""

import numpy as np
import sys
import os

paramfile = open("rad_mass_params.txt","w")

#algorithm begins:
#first use the corresponding rockstar bgc2 halo catalogue to find largest halo
data = np.loadtxt(sys.argv[1],skiprows=1) #halo catalogue data
blocks = sys.argv[2:len(sys.argv)]
largest_halo = np.zeros(6) #need mass, number of particles, x,y,z coordinates and virial radius (anything else?)
for i in range(0,len(data)): #scanning through all the halos
    if (data[i,3] > largest_halo[0]): #if it is the largest halo
        largest_halo[0] = data[i,3] #mass in Msun/h or Msun
        largest_halo[1] = data[i,4] #number of particles
        largest_halo[2] = data[i,11]*1000 #radius in pc/h
        largest_halo[3] = data[i,5]*1000 #x These in pc/h
        largest_halo[4] = data[i,6]*1000 #y
        largest_halo[5] = data[i,7]*1000 #z
        
#print to a parameter file the radius and particle mass
particlemass = largest_halo[0] / largest_halo[1] #total mass divided by number of particles, in Msun/h or Msun (figure it out)
paramfile.write(str(particlemass)+"\n")
paramfile.write(str(largest_halo[2])+"\n")


#now we have the location of the largest halo, so we can call our gadget block to lattice code to generate the .vis files
for i in blocks:
    print i
    command = "bsub -o "+i+".txt python readblock.py "+i+" " + str(largest_halo[3]) + " " + str(largest_halo[4]) + " " + str(largest_halo[5]) + " " + str(largest_halo[2])
    os.system(command)   

#NEXT WE NEED TO STORE THE LOCATIONS OF THE PARTICLES WITHIN THE RADIUS IN ANOTHER FILE, AND THEN SCAN THROUGH THOSE FOR CONCENTRIC RINGS

#build a 3d mesh again, but now the coordinates are scaled to be 10*r200 in each dimension from center of halo
#scan the particles in and place them in the corresponding bin, so only picking the ones within that square radii
#now we have the halo for visualization

#density profile:
#scan through all the halos within that radius and check location, mark which radii ring its in
#run back through and make the density profile, by summing up total in each concentric ring