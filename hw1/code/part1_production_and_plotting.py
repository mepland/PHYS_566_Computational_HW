#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import numpy as np
import matplotlib.pyplot as plt

# Produce the desired x ndarray, from -10.0 to 10.0 in steps of 0.05
x = np.linspace(-10.0, 10.0, 401)

# Declare the y ndarray
y1 = np.zeros(x.size)
y2 = np.zeros(x.size)

# Fill the y ndarrays
# Note that y2 at x = 0 is undefined not zero, so set it to None
for i in range(x.size):
	y1[i] = np.sin(x[i])
	if x[i] != 0.0:
		y2[i] = np.sin(x[i])/x[i]
	elif x[i] == 0.0:
		y2[i] = None

# Now plot and format
y1_plt, = plt.plot(x, y1, 'b', label='Sin(x)')
y2_plt, = plt.plot(x, y2, 'g', label='Sin(x)/x')

plt.axis([-11.0, 11.0, -1.1, 1.1])

plt.xlabel('x')
plt.ylabel('f(x)')

plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3, ncol=2, borderaxespad=0.0)


# Set output path, hard coded...
output_path = './output'

# Create the output dir, if it already exists don't crash, otherwise raise an exception
# Adapted from A-B-B's response to http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Note in python 3.4+ 'os.makedirs(output_path, exist_ok=True)' would handle all of this...
try: 
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise Exception('Problem creating output dir %s !!!\nA file with the same name probably already exists, please fix the conflict and run again.' % output_path)

# Print the plots
plt.savefig(output_path+'/part1.pdf')
plt.savefig(output_path+'/part1.eps')
plt.savefig(output_path+'/part1.png')
plt.savefig(output_path+'/part1.jpeg')

print 'Done!'

