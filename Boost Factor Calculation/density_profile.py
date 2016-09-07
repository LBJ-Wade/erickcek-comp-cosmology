# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 14:40:31 2015

@author: greensb
This script takes in as arguments the location of all particles in one file and the particle mass and halo radius from another file.

example:
python prof_comparison.py halo1_particles.txt halo1_params.txt
"""

import sys
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

def distance(x,y,z):
    return np.sqrt(x*x + y*y + z*z)

def nfw_log(r,rho_0,R_s,gamma):
    with np.errstate(invalid='ignore'): #SEE IF THIS IS STILL COMPLAINING
        rho = rho_0 / (((r / R_s)**gamma) * ((1 + (r / R_s))**(3.0 - gamma)))
    return np.log10(rho)
    
def nfw(r,rho_0,R_s,gamma):
    rho = rho_0 / (((r / R_s)**gamma) * ((1 + (r / R_s))**(3.0 - gamma)))
    return rho


res=128
paramfile = open(sys.argv[2],"rb") # contains particle mass and radius
partmass = np.float64(paramfile.readline())
comov_radius = np.float64(paramfile.readline()) #we are checking out to 10 x the radius (COMOVING)
sim_part_count = np.float64(paramfile.readline())
sim_size = np.float64(paramfile.readline())
redshift = 30.0#raw_input('Input final redshift: ') #figure out how to get this easily
scale = 1.0 / (1.0 + np.float64(redshift))

#read in our particles
particles = []
with open(sys.argv[1]) as particlefile:
particles = [line.split() for line in particlefile]
    

#here we make a density profile and get the generalized NFW profile from it for our integral
#below, we make a number of bins, and then run through the particles, finding their radius, and 
bin_count = res #number of bins (COULD USE REGULAR RES, OR DOUBLE RES, OR SOMETHING SO THAT INTEGRALS COMPUTED ON SAME RES...)
softening = (10**6)*.03*sim_size / sim_part_count#comov_radius/min_rad_frac #comoving softening length in pc/h
min_rad_frac = comov_radius / softening #greater than one, the number that comov_radius is divided by 
min_rad = comov_radius / min_rad_frac #minimum should be softening length
log_spaced_radii = np.logspace(np.log10(min_rad),np.log10(comov_radius),num=bin_count)
log_intervals = (np.log10(comov_radius) - np.log10(min_rad)) / (bin_count - 1)
bins = np.zeros(bin_count) #histogram essentially
for i in range(0,len(particles)):
x = np.float64(particles[i][0])
y = np.float64(particles[i][1])
z = np.float64(particles[i][2])
part_bin = np.ceil(np.log10(min_rad_frac * distance(x,y,z) / comov_radius) / log_intervals)
if(part_bin < 0):
    part_bin = 0 #just in case it is truly right on the edge of the max radius
#part_bin = np.floor((distance(x,y,z) / comov_radius) * bin_count)
if(part_bin >= bin_count):
    pass #we are working with a square, so we only want within radius data
#        elif (part_bin == bin_count):
#            part_bin = bin_count - 1
#            bins[part_bin] += 1
else:
    bins[part_bin] += 1

#print bins
#here we implement something tricky, which is that we collapse all inner bins such that the first bin has >100 particles in it
min_inner_count = 100
while(bins[0] < min_inner_count):
bins[1] += bins[0]
bins = np.delete(bins,0) #remove this index
log_spaced_radii = np.delete(log_spaced_radii,0)

bin_count = len(bins) #this is less than res if we have to collapse inner bins
radii = np.zeros(bin_count)
delta_rho = np.zeros(bin_count)
sigma = np.zeros(bin_count) 
true_density = np.zeros(bin_count)

for i in range(0,bin_count):
if(i == 0):
    delta_vol = (4.0/3.0 * np.pi)*((log_spaced_radii[0])**3)*(scale**3)
else:
    outer_vol = (4.0/3.0 * np.pi)*(log_spaced_radii[i]**3)*(scale**3)
    inner_vol = (4.0/3.0 * np.pi)*(log_spaced_radii[i-1]**3)*(scale**3)
    delta_vol = outer_vol - inner_vol
#outer_vol = (4.0/3.0 * np.pi)*(((i+1)* comov_radius / bin_count)**3)*(scale**3)
#inner_vol = (4.0/3.0 * np.pi)*((i* comov_radius / bin_count)**3)*(scale**3)
true_density[i] = bins[i] * partmass / delta_vol
delta_rho[i] = np.sqrt(bins[i]) * partmass / delta_vol
sigma[i] = 0.5*(np.abs(np.log10(true_density[i]+delta_rho[i])-np.log10(true_density[i]))+np.abs(np.log10(true_density[i]-delta_rho[i])-np.log10(true_density[i])))
#radii[i] = (i+1)*comov_radius*scale / bin_count
popt, pcov = opt.curve_fit(nfw_log, log_spaced_radii*scale, np.log10(true_density),None)
#rewrite better integrater over cartesian space
print sys.argv[1+(2*prof)].split('_')[0]
print "Optimized NFW profile parameters (rho_0,R_s,gamma):"
print popt
concentration = comov_radius*scale / popt[1] # r_vir/ r_S (radius is originally in comoving)
print "Concentration: " + str(concentration)
fit_density_plot = nfw_log(log_spaced_radii*scale,popt[0],popt[1],popt[2])
fit_density = nfw(log_spaced_radii*scale,popt[0],popt[1],popt[2])

#poisson errors sqrt(N particles in shell)/volume of the shell
plt.errorbar(np.log10(log_spaced_radii*scale),np.log10(true_density), yerr=sigma,label=sys.argv[1].split('_')[0])
plt.plot(np.log10(log_spaced_radii*scale),fit_density_plot)

plt.axvline(np.log10(softening*scale), label="Softening Length",color='green')
plt.legend(shadow = True)
plt.xlabel('Log10(Physical Radius) (pc/h)')
plt.ylabel('Density (M_sun/h)/(pc/h)^3')
plt.title('Density Profile')
plt.show()
