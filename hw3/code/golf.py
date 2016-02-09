#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

########################################################
# Set fixed parameters, SI Units
mass = 46*0.001 # mass of ball = 46 grams
v0 = 70.0 # Initial velocity
rho = 1.29 # Density of air (at sea level)
A =  0.0014 # Frontal Area
g = 9.80665 # g at sea level

theta0 = [45, 30, 15, 9] # Initial launch angles

C = [0.5, 7.0] # Drag Coefficient TODO units are ??

dt = 0.0005 # Time step of simulation, empirically set to satisfy max_r
max_r = 0.05 # Maximum distance allowed between time steps
# For comparison, Solving 0.0014 = 2 pi r^2 for half the surface area of a sphere gives us the balls radius r ~= 0.015


########################################################
# Print out starting values
print '\nMass is: %.3f' % mass
print 'Initial Velocity is: %.1f' % v0
print 'Air Density Rho is: %.2f' % rho
print 'Frontal Area A is: %.4f' % A
print 'g is: %.5f' % g

print '\nLaunch Angles Theta are: '
print theta0

print '\nDrag Coefficients are: '
print C

print '\nSimulation Time Step is: %.2f' % dt
print 'Simulation Max Step Size is: %.2f' % max_r

########################################################
# Define time_step class to hold all the relevant parameters at a time step
class time_step:

    def __init__(self, time):
	self.t = time
        self.x = -99.0
	self.y = -99.0
	self.vx = -99.0
	self.vy = -99.0
	
	# Return theta of the time step when called
    def theta(self):
	if self.x == -99.0 and self.y == -99.0 and self.vx == -99.0 and self.vy == -99.0: print 'WARNING COMPUTING THETA ON DEFAULT TIME STEP'
	return np.arctan2(self.vy, self.vx)

	# Return spatial separation r with another times tep
    def r(self, time_step2):
	if self.x == -99.0 and self.y == -99.0 and self.vx == -99.0 and self.vy == -99.0: print 'WARNING COMPUTING r WITH ON DEFAULT TIME STEP'
	return np.sqrt( np.square(self.x - time_step2.x) + np.square(self.y - time_step2.y))


########################################################
# TODO Define a function to run the simulation based on which case we are in and what theta0 we want



# Start the list of time_steps
run = [time_step(0.0)]

# Set the initial values
run[0].x = 0.0
run[0].y = 0.0
run[0].vx = v0*np.cos(theta0[0]) # TODO change theta0 in overarching loop
run[0].vy = v0*np.sin(theta0[0])

# Loop until we hit the ground
current_y = 99.0 # use to halt while loop
i = 0 # current time step

while current_y > 0.0:
	i = i + 1 # increment time step
	run.append(time_step(i*dt)) # append a new time_step object

	# Compute the new position
	run[i].x = run[i-1].x + run[i-1].vx*dt
	run[i].y = run[i-1].y + run[i-1].vy*dt
	# Compute the new velocity TODO change forces here
	run[i].vx = run[i-1].vx
	run[i].vy = run[i-1].vy + (dt/mass)*(-mass*g)

	# Make sure we didn't move to far, notify user and exit if we did
	if(run[i].r(run[i-1]) > max_r):
		print 'On time step %d the ball moved to far, r = %.3f' % (i, run[i].r(run[i-1]))
		sys.exit()
	current_y = run[i].y # update current_y

########################################################
# Create ndarrays that we can plot

# Declare the x and y ndarrays based on our finished run's size
x_ndarray = np.zeros(len(run))
y_ndarray = np.zeros(len(run))

# Fill the ndarrays
for j in range(len(run)):
	x_ndarray[j] = run[j].x
	y_ndarray[j] = run[j].y


########################################################
# Create the plot object
m_label = '$\\theta_{0}$ = %.1f' % theta0[0] # TODO alter in big loop

plt.scatter(x_ndarray, y_ndarray, marker='.', label=m_label)




########################################################
# Format the plot

# Warning, y axis range is hard coded here!
# plt.axis([0, stop_time, 3*10**8, 6*10**9])

plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')

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



print '\n\nDone!\n'

