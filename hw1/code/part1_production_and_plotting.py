#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Produce the desired x vector, from -10.0 to 10.0 in steps of 0.05
x = np.linspace(-10.0,10.0,401)

# Declare the y vectors
y1 = np.zeros(x.size)
y2 = np.zeros(x.size)

# Fill the y vectors
# Note that y2 at x = 0 is undefined not zero, but will be zero in this scheme TODO
for i in range(x.size):
	y1[i] = np.sin(x[i])
	if x[i] != 0.0:
		y2[i] = np.sin(x[i])/x[i]
	elif x[i] == 0.0:
		y2[i] = 0.0 # Already should be 0.0 from initialization, but may change this later

# Now plot x, y1, y2 all together
y1_plt, = plt.plot(x, y1, 'b', label='Sin(x)')
y2_plt, = plt.plot(x, y2, 'g', label='Sin(x)/x')

plt.axis([-11, 11, -1.1, 1.1])

plt.xlabel('x')
plt.ylabel('f(x)')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, borderaxespad=0.)

plt.savefig('part1.pdf')
plt.savefig('part1.eps')
plt.savefig('part1.png')
plt.savefig('part1.jpeg')

#print "\n y.size returns %i" % y.size

# write output to file
#output = open("output.txt","w")
#for i in range(y.size):
#	output.write('%2.6f \n' % y[i])

# close output file
#output.close()




