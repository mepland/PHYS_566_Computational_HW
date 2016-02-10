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

vartheta0 = [45, 30, 15, 9, 60] # Initial launch angles (degrees)
# Convert to radians
for k in range(len(vartheta0)):
	vartheta0[k] = vartheta0[k]*(np.pi/180.0)

# Drag Coefficients and transition velocity
C = 0.5
C_prime = 7.0
v_transition = 14.0

dt = 0.0005 # Time step of simulation, empirically set to satisfy max_r
max_r = 0.05 # Maximum distance allowed between time steps
# For comparison, Solving 0.0014 = 2 pi r^2 for half the surface area of a sphere gives us the balls radius r ~= 0.015


########################################################
# Print out starting values
print '\nBeginning golf.py simulation'
print 'All units are SI unless otherwise specified'

print '\nMass is: %.3f' % mass
print 'Initial Velocity is: %.1f' % v0
print 'Air Density Rho is: %.2f' % rho
print 'Frontal Area A is: %.4f' % A
print 'g is: %.5f' % g

print '\nLaunch Angles Vartheta are: '
for k in range(len(vartheta0)):
	print '%.4f radians, %.1f degrees' % (vartheta0[k], vartheta0[k]*(180.0/np.pi))

print '\nDrag Coefficients are C = %.1f and C\' = %.1f' % (C, C_prime)
print 'Transition Velocity is: %.1f' % v_transition

print '\nSimulation Time Step is: %.5f' % dt
print 'Simulation Max Step Size is: %.3f' % max_r

# TODO add ---s

########################################################
# Define time_step class to hold all the relevant parameters at a time step
class time_step:

    def __init__(self, time):
	self.t = time
        self.x = -99.0
	self.y = -99.0
	self.vx = -99.0
	self.vy = -99.0
	
	# Return vartheta of the time step when called
    def vartheta(self):
	if self.x == -99.0 and self.y == -99.0 and self.vx == -99.0 and self.vy == -99.0: print 'WARNING COMPUTING THETA ON DEFAULT TIME STEP'
	return np.arctan2(self.vy, self.vx)

	# Return the speed v of the time step when called
    def vmag(self):
	if self.x == -99.0 and self.y == -99.0 and self.vx == -99.0 and self.vy == -99.0: print 'WARNING COMPUTING VMAG ON DEFAULT TIME STEP'
	return np.sqrt(np.square(self.vx) + np.square(self.vy))

	# Return spatial separation r from another time step
    def r(self, time_step2):
	if self.x == -99.0 and self.y == -99.0 and self.vx == -99.0 and self.vy == -99.0: print 'WARNING COMPUTING r WITH ON DEFAULT TIME STEP'
	return np.sqrt( np.square(self.x - time_step2.x) + np.square(self.y - time_step2.y))
# end class for time_step

########################################################
# Define a function to run the simulation based on which case we are in and what vartheta0 we want
def run_sim(case, m_vartheta0, m_color):
	
	# Start the list of time_steps
	run = [time_step(0.0)]
	
	# Set the initial values
	run[0].x = 0.0
	run[0].y = 0.0
	run[0].vx = v0*np.cos(m_vartheta0)
	run[0].vy = v0*np.sin(m_vartheta0)

	# Loop until we hit the ground
	current_y = 99.0 # use to halt while loop
	i = 0 # current time step
	C_dimpled = C # Declare C_dimpled for use later, just set it to C for now

	while current_y > 0.0:
		i = i + 1 # increment time step
		run.append(time_step(i*dt)) # append a new time_step object
	
		# Compute the new position
		run[i].x = run[i-1].x + run[i-1].vx*dt
		run[i].y = run[i-1].y + run[i-1].vy*dt
		# Compute the new velocity, see writeup for physics...
		if(case == 'a'): # ideal
			run[i].vx = run[i-1].vx + (dt/mass)*( 0.0 )
			run[i].vy = run[i-1].vy + (dt/mass)*( -mass*g )
		if(case == 'b'): # smooth ball with drag
			F_drag_smooth = -C*rho*A*np.square(run[i-1].vmag())
			run[i].vx = run[i-1].vx + (dt/mass)*( F_drag_smooth*np.cos(run[i-1].vartheta()) )
			run[i].vy = run[i-1].vy + (dt/mass)*( F_drag_smooth*np.sin(run[i-1].vartheta()) -mass*g )
		if(case == 'c'): # dimpled ball with drag
			# see what v regime we are in and pick C accordingly
			if( run[i-1].vmag() <= v_transition):
				C_dimpled = C # 1/2
				# print 'v = %.2f, C_dimpled = %.3f, BELOW TRANSITION' % (run[i-1].vmag(), C_dimpled)
			else:
				C_dimpled = C_prime/run[i-1].vmag() # C'/v
				# print 'v = %.2f, C_dimpled = %.3f' % (run[i-1].vmag(), C_dimpled)
			F_drag_dimpled = -C_dimpled*rho*A*np.square(run[i-1].vmag())
			run[i].vx = run[i-1].vx + (dt/mass)*( F_drag_dimpled*np.cos(run[i-1].vartheta()) )
			run[i].vy = run[i-1].vy + (dt/mass)*( F_drag_dimpled*np.sin(run[i-1].vartheta()) -mass*g )
		# TODO add more cases
		# Make sure we didn't move to far, notify user and exit if we did
		if(run[i].r(run[i-1]) > max_r):
			print 'On time step %d the ball moved to far, r = %.3f\nProgram Exiting' % (i, run[i].r(run[i-1]))
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
	m_label = '$\\vartheta_{0}$ = %.1f$^{\circ}$' % (m_vartheta0*(180.0/np.pi))
	print 'Run Completed'
	return plt.scatter(x_ndarray, y_ndarray, marker='.', label=m_label, c=m_color, edgecolors='none')
# end def for run_sim

########################################################
# Define a function to run the simulation multiple times for all the vartheta0's 
def loop_vartheta0(case):
	# Define the colors we want to use, include some extra to be safe
	colors = ['black', 'blue', 'green', 'cyan', 'magenta', 'crimson']
	scatters = []
	for k in range(len(vartheta0)):
		print 'Running with vartheta0 = %.1f degrees' % (vartheta0[k]*(180.0/np.pi))
		scatters.append(run_sim(case, vartheta0[k], colors[k]))
	return scatters
# end def for loop_vartheta0

########################################################
# Define a function to run, plot, and print a whole case 
def run_case(case, title, fname):
	print '\nRunning Case: %s' % title
	print '-----------------------------------------\n'

	fig = plt.figure(case) # get a separate figure for the case
	ax = fig.add_subplot(111) # Get the axes, effectively

	loop_vartheta0(case) # generate the scatters

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
	# fig.savefig(output_path+'/'+fname+'.pdf')
	# Due to the number of points pdfs will be VERY large
	# So instead only make high quality pngs
	fig.savefig(output_path+'/'+fname+'.png', dpi=900)

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
#'''
run_case('a', 'Ideal Ball', 'ideal')
run_case('b', 'Smooth Ball with Drag', 'smooth_ball_w_drag')
run_case('c', 'Dimpled Ball with Drag', 'dimpled_ball_w_drag')
#'''
#run_case('d', 'Dimpled Ball with Drag and Spin', 'dimpled_ball_w_drag_and_spin')


########################################################
print '\n\nDone!\n'

