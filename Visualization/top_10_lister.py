#takes as an argument a halo catalog and then prints out the information about
#the top 10 largest in the catalog

import numpy as np
import sys

data=np.loadtxt(sys.argv[1])
data = data[data[:,1].argsort()]
data = data[len(data)-10:,:]

print data
