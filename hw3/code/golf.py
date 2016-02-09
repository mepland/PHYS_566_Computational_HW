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

theta0 = [45, 30, 15, 9, 60] # Initial launch angles (degrees)
# Convert to radians
for k in range(len(theta0)):
	theta0[k] = theta0[k]*(np.pi/180.0)

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
for k in range(len(theta0)):
	print '%.4f radians, %.1f degrees' % (theta0[k], theta0[k]*(180.0/np.pi))

print '\nDrag Coefficients are: '
print C

print '\nSimulation Time Step is: %.5f' % dt
print 'Simulation Max Step Size is: %.3f\n' % max_r

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
# end class for time_step

########################################################
# Define a function to run the simulation based on which case we are in and what theta0 we want
def run_sim(case, m_theta0, m_color):
	
	# Start the list of time_steps
	run = [time_step(0.0)]
	
	# Set the initial values
	run[0].x = 0.0
	run[0].y = 0.0
	run[0].vx = v0*np.cos(m_theta0)
	run[0].vy = v0*np.sin(m_theta0)

	# Loop until we hit the ground
	current_y = 99.0 # use to halt while loop
	i = 0 # current time step

	while current_y > 0.0:
		i = i + 1 # increment time step
		run.append(time_step(i*dt)) # append a new time_step object
	
		# Compute the new position
		run[i].x = run[i-1].x + run[i-1].vx*dt
		run[i].y = run[i-1].y + run[i-1].vy*dt
		# Compute the new velocity
		if(case == 'a'):
			run[i].vx = run[i-1].vx
			run[i].vy = run[i-1].vy + (dt/mass)*(-mass*g)
		#if(case == 'b'):# TODO
			#run[i].vx = run[i-1].vx
			#run[i].vy = run[i-1].vy + (dt/mass)*(-mass*g)

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
	# Create and return the scatter object
	m_label = '$\\theta_{0}$ = %.1f$^{\circ}$' % (m_theta0*(180.0/np.pi))
	print 'Run Completed'
	return plt.scatter(x_ndarray, y_ndarray, marker='.', label=m_label, c=m_color, edgecolors='none')
# end def for run_sim

########################################################
# Define a function to run the simulation multiple times for all the theta0's 
def loop_theta0(case):
	# Define the colors we want to use, include some extra to be safe
	colors = ['black', 'blue', 'green', 'cyan', 'magenta', 'crimson']
	scatters = []
	for k in range(len(theta0)):
		print 'Running with theta0 = %.1f degrees' % (theta0[k]*(180.0/np.pi))
		scatters.append(run_sim(case, theta0[k], colors[k]))
	return scatters
# end def for loop_theta0

########################################################
# Define a function to run, plot, and print a whole case 
def run_case(case, title, fname):
	print 'Running Case: %s' % title
	print '-----------------------------------\n'

	fig = plt.figure(case) # get a separate figure for the case
	ax = fig.add_subplot(111) # Get the axes, effectively

	loop_theta0(case) # generate the scatters

	##########################
	# Format the plot

	ax.set_title(title)

	ax.set_xlabel('$x$ [m]')
	ax.set_ylabel('$y$ [m]')

	ax.set_xlim((0.0, None))
	ax.set_ylim((0.0, None))

	plt.legend()

	##########################
	# Print the plot
	fig.savefig(output_path+'/'+fname+'.pdf')
	fig.savefig(output_path+'/'+fname+'.png')

	print 'Case Completed'
# end def for run_case



########################################################
# Create the output dir, if it already exists don't crash, otherwise raise an exception
# Adapted from A-B-B's response to http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Note in python 3.4+ 'os.makedirs(output_path, exist_ok=True)' would handle all of this...
output_path = './output' # Set output path, hard coded...
try: 
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise Exception('Problem creating output dir %s !!!\nA file with the same name probably already exists, please fix the conflict and run again.' % output_path)

########################################################
# Finally, actually run over the four cases

run_case('a', 'Ideal', 'ideal')
# run_case('b', 'Smooth Ball with Drag', 'smooth_ball_w_drag')




########################################################
print '\n\nDone!\n'

