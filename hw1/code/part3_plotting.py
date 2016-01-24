#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import numpy as np
import matplotlib.pyplot as plt

# Set input path, hard coded...
input_path = './output'

# open input file and see how many data points/lines there are
input_file = open(input_path+'/data.txt','r')
num_data = sum(1 for line in input_file)
#print 'Input file %s has %i lines' % (input_file.name, num_data)

# Declare the y ndarray
x  = np.zeros(num_data)
y1 = np.zeros(num_data)
y2 = np.zeros(num_data)

# Return to the top of the file and read the data in
input_file.seek(0)

for i in range(num_data):
	line = input_file.readline()
	x[i]  = line.split('\t')[0]
	y1[i] = line.split('\t')[1]
	y2[i] = line.split('\t')[2]
	#print '%(x).2f \t %(y1).6f \t %(y2).6f' % {'x': x[i], 'y1': y1[i], 'y2': y2[i]}

# Close the input file
input_file.close()


#######################################################
# Recycle the plotting code from part 1

# Now plot and format
y1_plt, = plt.plot(x, y1, 'b', label='Sin(x)')
y2_plt, = plt.plot(x, y2, 'g', label='Sin(x)/x')

plt.axis([-11.0, 11.0, -1.1, 1.1])

plt.xlabel('x')
plt.ylabel('f(x)')

plt.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3, ncol=2, borderaxespad=0.0)


# Set output path, hard coded...
output_path = './output'

# Create the output dir, if it exists don't crash, otherwise raise an exception
try: 
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise Exception('Problem creating output dir %s !!!\nA file with the same name probably already exists, please fix the conflict and run again.' % output_path)

# Print the plots
plt.savefig(output_path+'/part3.pdf')
plt.savefig(output_path+'/part3.eps')
plt.savefig(output_path+'/part3.png')
plt.savefig(output_path+'/part3.jpeg')

print 'Done!'

