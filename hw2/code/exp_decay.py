#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import numpy as np
import matplotlib.pyplot as plt


# Set inital parameters
# Units: Time ~ years, N ~ raw number, mass ~ kg
half_life = 5700
initial_mass = 10**-12
#stop_time = 20000
stop_time = 30

# Compute initial N0 from initial mass
NA = 6.022*(10**23)
N0 = NA*((initial_mass*1000)/14)

# Print out starting values
print 'Half-Live is: %d' % half_life
print 'Initial Mass is: %2.2e' % initial_mass
print 'Initial Number (N0) is: %2.2e' % N0
print 'Stop Time is: %d' % stop_time

'''
# Debugging TODO Delete when submitting
print type(half_life)
print type(initial_mass)
print type(NA)
print type(N0)
'''

# Set time steps
dt1 = 10
dt2 = 100
dt3 = 1000

# Produce the desired t ndarrays, from 0.0 to stop_time in steps of dt
# need to add stop_time+dt to get stop_time
t1 = np.arange(0,stop_time+dt1,dt1)

# Declare the N ndarrays
N1 = np.zeros(t1.size)




'''
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
plt.savefig(output_path+'/part1.png')
'''
print 'Done!'

