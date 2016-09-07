# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 23:53:39 2015

@author: greensb
"""

import numpy as np
import sys

#first use the corresponding rockstar bgc2 halo catalogue to find largest halo

data = np.loadtxt(sys.argv[1]) #halo catalogue data
top100_subhalos = np.zeros(100,dtype=np.int) #we need as many as there are 
top100 = np.loadtxt(sys.argv[2]) #list of top100 halos
for i in range(0,len(data)): #scanning through all the halos
    parent = data[i,12]
    while(parent != -1):
        if(parent in top100): #if the parent is our halo in question
            index = np.where(top100 == parent)[0]
            top100_subhalos[index] += 1
        parent_index = np.where(data[:,0] == parent)[0]
        parent = data[parent_index,12]
        
avg = np.mean(top100_subhalos)
std = np.std(top100_subhalos)


print avg
print std
print "max subhalos"
print max(top100_subhalos)