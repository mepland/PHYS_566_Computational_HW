#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python

#from pylab import *
import numpy as np
#from scipy import *


# Produce the desired vector
x = np.linspace(0,10,101)

print x

y = np.sin(x)

print y

#print "\n y.size returns %i" % y.size

# write output to file
output = open("output.txt","w")
for i in range(y.size):
	output.write('%2.6f \n' % y[i])


# close output file
output.close()




