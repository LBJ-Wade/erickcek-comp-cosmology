import numpy as np
import matplotlib.pyplot as plt
import sys

#Sets the minimum width of each bin and also the amount by width a bin grows if necessary.
#These values can be changed to tweak the output. As a general guidline, it seems to work well if delta is about 1/4 of min width
#Decreasing these values -> smaller bins, while increasing -> larger bins
#Three different sized sets have been provided
minBinWidth = .4
deltaBinWidth = 0.1
#minBinWidth = .2
#deltaBinWidth = 0.05
#minBinWidth = 0.8
#deltaBinWidth = 0.2

##############################################################################################################
#USER INPUT SECTION#
pltType = 'scat'

#User must enter the number of files to be processed
numFiles = len(sys.argv)-1
#modify this to plot for comparison

pltTime = [0] * numFiles

#Asks user for the size of the box. Must be given to compare mass functions across different box sizes
boxSize = float(raw_input('''Enter the box size used in this simulation in Mpc:  '''))

#The value of h for the simulation. Usually around 0.72
hval = float(raw_input('''Enter the value of h used in this simulation:  '''))

#eventually can ask user the amount of particles and then determine the PARTICLE MASS for minimum size VERTICAL LINE

#Divides box size by h to scale appropriately 
boxSize = boxSize/hval
#Volume of the box is boxSize cubed
vol = boxSize**3


######################################################################################################################################################################
#MAIN CODE BODY#

#This is the main for statement of the code, which tells us that the code will loop numFiles times
for j in xrange (1, numFiles+1):
  #Load the data file. All files must end in .ascii at the moment, although this could easily be changed by asking the user what kind of file they are inputting.
  data = np.loadtxt(sys.argv[j])
  #remove subhalos
  data = data[data[:,12] == -1] #only host halos
  #Rockstar gives .ascii outputs with the list of virial masses in the third column. If you would like to use a different column, change the line below
  masses = data[:,2] / hval
  #This is part of the dynamic binning routine. Based on the number of halos found, this sets the minimum number of halos that must be found in each bin.
  #Some quick experimentation showed that this relationship is logarithmic. The current formulation (10.54 * np.log(len(masses)) - 40.36) is based on a small number of points
  #Feel free to further refine this formula, it may well give better results!
  #NOTE Currently, this does not work well for small resolution (N=32), but does work well for higher res
  minHalosInBin = 10.54 * np.log(len(masses)) - 40.36
  #The code now logs and sorts the masses. This is done for ease of manipulation and mass function determination.
  masses = np.log(masses)
  masses = np.sort(masses)
  #Because the array has been sorted, the min and max are just the first and last elements.
  minMass = masses[0]
  maxMass = masses[-1]
      
  #This begins the core routine. currentLoc describes where the code thinks that the next right hand side of a bin will be. It's initial value is slightly strange looking.
  #minMass represents the lowest mass value, but we don't want to start right on that value. So we take that value, and subtract the small delta value that we define above. 
  #That value will actually be the start point, the left most edge of the first bin. You can see this as it is the first value in binSeq, the list/array that will eventually define
  #All the edges of all of the bins. Because we want currentLoc to be the right hand edge, we add the minBinWidth to it. So we now have a left edge slightly lower than the lowest mass
  #value and a right edge that is exactly 1 minBinWidth away. Great!
  currentLoc = minMass-deltaBinWidth+minBinWidth
  binSeq = [currentLoc-minBinWidth]
  i=0
  binCount =0
  #The first while checks to make sure that our current expected right edge is not beyond where we want the rightmost edge of the last bin to possibly be. Again, we set this to slightly
  #higher than the largest value. Mostly this is just done for safety.
  while currentLoc < maxMass + deltaBinWidth:  
    #You run into a problem if the last bin never reaches the min bin count, so this just checks to see if the current location is near enough to the final bin edge to just
    #lump the rest of the halos into the final bin. This value seems to work pretty well.
    if (maxMass+deltaBinWidth) - currentLoc < .1:
    #At that point, we make currentLoc the final bin edge and then append that edge to the list. If you ever notice that the last bin is particularly large, this is probably why.
      currentLoc = maxMass + deltaBinWidth
      binSeq.append(currentLoc)
    else:
    #If we are not close to the end, then we want to count through the list of masses that are contained in this bin.
      while masses[i]<currentLoc:
      #So if the current mass is lower than the expected right wall, add 1 to the bin count and the overall counter. Once it has counted all the halos in the range, move on.
        binCount+=1
        i+=1    
      if binCount >= minHalosInBin: #Checks to see if the number of halos counted in the bin meets the minimum halos requirement
	binCount = 0 #If it does, reset the bin count
	binSeq.append(currentLoc) #Add the currentLocation to the list of bin edges
	currentLoc += minBinWidth #And then increase the current location by the min width, creating a new prospective bin
      else: #If we don't have enough halos, increase the width of the bin by a small number and count how many halos are newly included!
        currentLoc += deltaBinWidth
	      
  numBins = len(binSeq)-1 #The number of bins is dynamic, so we need to calculate what it is
  binSeq = np.array(binSeq) #convert binSeq to an array to make it easy to use
  binSize = np.array([0.0]*numBins) # We need bin size for the mass function, so we calculate it here. It is pretty simple, just takes the difference in the edges.
  for x in xrange(0, numBins):
    binSize[x] = binSeq[x+1] - binSeq[x]
  #Numpys histogram command documentation can be found online, but here it outputs the binEdgs (redundant, because we basically fed it the bin Edges), and the halosInBin (not redundant, becuase we never actually counted them above
  halosInBin,binEdges=np.histogram(masses,binSeq)
  binCenters = 0.5*(binSeq[1:]+binSeq[:-1])
  poissonError = np.sqrt(halosInBin) #Poisson error, taking sqrt of the halos in each bin
  binCenters = np.log10(np.exp(binCenters)) #Leading edges have been logged, so we unlog and then take log10 because we want log10 axes
  massFunc = np.log10(np.divide(halosInBin, (binSize*vol))) #Mass function is N/(vol*binSize) which gives units of Mpc^-3
  #Below we handle assymetric error calculations. We want to calculate error relative to the data, so we subtract the mass function from the calculation of error positions.
  #This seems complicated and redundant, but is necessary because of the log scale. Recall that log(x-y) is not the same as log(x) - log(y)
  upperError = np.log10(np.divide(halosInBin, binSize*vol) + np.divide(poissonError, binSize*vol)) - massFunc
  lowerError = np.log10(np.divide(halosInBin, binSize*vol) - np.divide(poissonError, binSize*vol)) - massFunc
  #Take absolute value of the errors. This is done because of how pyplot's error bars code works. Explained below, or by pyplots documentation online
  upperError = np.abs(upperError)
  lowerError = np.abs(lowerError)
  #We make a 1x2 array of 1d arrays. Pyplot takes these values and assumes that the first array is your low error bound, so it plots -column1 and then +column 2 relative to the error.
  #This is why we must take the absolute value above. Lower error comes out negative, but pyplot assumes it will be a positive value, so we need to fix that.
  assymError = [lowerError, upperError]
  
  #Below is the plot command control for the three plot types. Plt label will eventually be used to label the different lines on the graph.
  if pltType == 'scat':
    plt.figure(j) #If the user wants to plot multiple histograms, they need to be separate figures. This handles that, creating a figure with the index of the for loop above.
    pltLabel = sys.argv[j].split('.')[0]

  #The actual plot command. Plots massFunc vs binCenters (center of the mass bins, serves as x location), the y axis error, and the labels.
  #the fmt and ecolor simply define that it should plot connected circles with red error bars
  plt.errorbar(binCenters, massFunc, yerr = assymError, fmt='-o', ecolor = 'r', label = pltLabel)
  plt.legend(shadow = True)#Draws the legend
  #Axis labels
  plt.xlabel('Log10(Mass)') #check the units!!!!
  plt.ylabel('log10(dN/dlnM)')
  plt.draw()#Draws the figure, but does not show it. That command is held until the end.
  
  
  
#bring in our comparison data from PS
data_compare = np.loadtxt("massfunction-PS-pl15-6000-t30mc40.dat")
t = np.log10(data_compare[:,0]/hval)
s = np.log10(data_compare[:,1])
plt.plot(t,s, label="Press-Schecter")

#bring in our comparison data from ST
data_compare = np.loadtxt("massfunction-ST-pl15-6000-t30mc40.dat")
t1 = np.log10(data_compare[:,0]/hval)
s1 = np.log10(data_compare[:,1])
plt.plot(t1,s1,label="Sheth-Tormen")

part_mass = 1.72455426627e-11 /hval #hardcoded for .000030, 512 runs
log_min_mass = np.log10(part_mass*100)
plt.axvline(log_min_mass, label="Min. Trust Mass")
plt.legend(shadow = True)#Draws the legend

plt.legend(shadow = True)#Draws the legend
plt.draw()#Draws the figure, but does not show it. That command is held until the end.

plt.show()#Shows the final plot(s), with multiple lines overplotted if chosen
