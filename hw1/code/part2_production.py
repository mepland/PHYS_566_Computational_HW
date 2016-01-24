#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import numpy as np
import matplotlib.pyplot as plt

# Produce the desired x ndarray, from -10.0 to 10.0 in steps of 0.05
x = np.linspace(-10.0, 10.0, 401)

# Declare the y ndarrays
y1 = np.zeros(x.size)
y2 = np.zeros(x.size)

# Fill the y ndarrays
# Note that y2 at x = 0 is undefined not zero, so set it to None TODO watch for crashes in part 3 on nan
for i in range(x.size):
	y1[i] = np.sin(x[i])
	if x[i] != 0.0:
		y2[i] = np.sin(x[i])/x[i]
	elif x[i] == 0.0:
		y2[i] = None

# Set output path, hard coded...
output_path = './output'

# Create the output dir, if it exists don't crash, otherwise raise an exception
try:
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise Exception('Problem creating output dir %s !!!\nA file with the same name probably already exists, please fix the conflict and run again.' % output_path)

# write output to file
output_file = open(output_path+'/data.txt','w')
for i in range(x.size):
        output_file.write('%(x).2f \t %(y1).6f \t %(y2).6f \n' % {'x': x[i], 'y1': y1[i], 'y2': y2[i]})

# close output file
output_file.close()

print 'Done!'

