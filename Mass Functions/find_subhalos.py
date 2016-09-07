# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 23:53:39 2015

@author: greensb
"""

import numpy as np
import sys

#first use the corresponding rockstar bgc2 halo catalogue to find largest halo


data = np.loadtxt(sys.argv[1]) #halo catalogue data
halo = int(sys.argv[2]) #halo id of host halo
sub_halos = 0
for i in range(0,len(data)): #scanning through all the halos
    parent = data[i,12]
    while(parent != -1):
        if(parent == halo): #if the parent is our halo in question
            sub_halos += 1
        parent_index = np.where(data[:,0] == parent)[0]
        parent = data[parent_index,12]
        
        
print sub_halos