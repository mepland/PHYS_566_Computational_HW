#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import numpy as np
import matplotlib.pyplot as plt

########################################################
# Set inital parameters
# Units: Time ~ years, N ~ raw number, mass ~ kg
half_life = 5700.0
initial_mass = 10**-12
stop_time = 20000

########################################################
# Compute initial N0 from initial mass
NA = 6.022*(10**23)
N0 = NA*((initial_mass*1000.0)/14.0)

# Compute decay constant tau from half-life
tau = half_life/np.log(2)

# Print out starting values
print '\nHalf-Live is: %d' % half_life
print 'Decay Constant is: %d' % tau
print 'Initial Mass is: %2.2e' % initial_mass
print 'Initial Number (N0) is: %2.2e' % N0
print 'Stop Time is: %d' % stop_time

########################################################
# Set time steps
dt1 = 10.0
dt2 = 100.0
dt3 = 1000.0

# Produce the desired t ndarrays, from 0.0 to stop_time in steps of dt
# need to add stop_time+dt to get stop_time
t1 = np.arange(0,stop_time+dt1,dt1)
t2 = np.arange(0,stop_time+dt2,dt2)
t3 = np.arange(0,stop_time+dt3,dt3)

# Declare the N ndarrays
N1 = np.zeros(t1.size)
N2 = np.zeros(t2.size)
N3 = np.zeros(t3.size)

# Initialize the N ndarrays to N0
N1[0] = N0
N2[0] = N0
N3[0] = N0

# Fill the ndarrays
for i in range(1,t1.size):
	N1[i] = (1 - (dt1/tau))*N1[i-1]

for i in range(1,t2.size):
	N2[i] = (1 - (dt2/tau))*N2[i-1]

for i in range(1,t3.size):
	N3[i] = (1 - (dt3/tau))*N3[i-1]

########################################################
# Create the plots
R1_label = '$\Delta t = %.0d$' % dt1
R1_plt, = plt.plot(t1, N1/tau, 'b--', label=R1_label)

R2_label = '$\Delta t = %.0d$' % dt2
R2_plt, = plt.plot(t2, N2/tau, 'g-.', label=R2_label)

R3_label = '$\Delta t = %.0d$' % dt3
R3_plt, = plt.plot(t3, N3/tau, 'm:', label=R3_label)


R_explicit_plt, = plt.plot(t1, (N0/tau)*np.exp(-t1/tau), 'k', label='$R\left(t\\right) = \left(N_{0}/\\tau\\right) e^{\left(-t/\\tau\\right)}$')
m_dashes = [4, 8, 8, 8]  # number of points: on, off, on, off
R_explicit_plt.set_dashes(m_dashes)


########################################################
# Format the plot

# Warning, y axis range is hard coded here!
plt.axis([0, stop_time, 3*10**8, 6*10**9])

plt.xlabel('$t$ [Years]')
plt.ylabel('$R\left(t\\right)$ [Counts / Year]')

plt.legend()

########################################################
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

########################################################
# Print the plot
plt.savefig(output_path+'/plot.pdf')
plt.savefig(output_path+'/plot.png')

########################################################
# Make y axis log scale
plt.yscale('log')

# Re-Print the plot
plt.savefig(output_path+'/log_plot.pdf')
plt.savefig(output_path+'/log_plot.png')

########################################################
# Create diagnostic plot to check the half-life is in the correct place...

# Add vertical line at T1/2
plt.axvline(x=half_life, color='r', linestyle='-', label='$T_{1/2}$')

# Add horizontal line at (N0/2)/tau
plt.axhline(y=(N0/2)/tau, color='c', linestyle='-', label='$N_{0}/(2$$\\tau)$')

# Redraw the legend
plt.legend()

# Print the diagnostic plot
plt.savefig(output_path+'/diagnostic_plot.pdf')
plt.savefig(output_path+'/diagnostic_plot.png')


########################################################
# TODO Compute and printout the desirec accuracy calculations 



print '\nDone!\n'

