#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl # We need to import all of matplotlib just to set rcParams once...
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

########################################################
# Set fixed/global parameters [SI units]

Q_over_epsilon0 = 1.0 # point charge Q over epsilon0 [Vm]

# TODO Run values
a = 0.6 # dipole separation length [m]
R = 10.0 # Location of Spherical Boundary Condition [m]

# TODO debuging values
#a = 2.0
#R = 3.0

########################################################
# Print out fixed values
print '\nBeginning poisson_dipole.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nQ over epsilon0 = %.1f [Vm]' % Q_over_epsilon0
print 'dipole seperation a = %.1f [m]' % a
print 'Circular boundary R = %.1f [m]' % R

print '\n---------------------------------------------'
print '---------------------------------------------\n'

########################################################
########################################################

########################################################
# Define a function to create the output dir
# If it already exists don't crash, otherwise raise an exception
# Adapted from A-B-B's response to http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Note in python 3.4+ 'os.makedirs(output_path, exist_ok=True)' would handle all of this...
def make_path(path):
	try: 
	    os.makedirs(path)
	except OSError:
	    if not os.path.isdir(path):
	        raise Exception('Problem creating output dir %s !!!\nA file with the same name probably already exists, please fix the conflict and run again.' % output_path)
# end def for make_path

########################################################
# Define a function to run the simulation
def run_sim(L, Dx, sim_method, halt_method, epsilon, fixed_accuracy, debugging):
	print 'Beginning run_sim:'	

	# Find the index of the center/origin
	l_center = (L-1)/2 # +1 to get to center, -1 as we start from 0
	
	########################################################
	# Define a function to take ndarray indexes and return spatial coordinates of the BIN CENTERS 
	# TODO move out of this function?
	def find_spatial(y_index,x_index):
		# X: 0 ~ - R, l_center=(L-1)/2 ~ 0.0, L-1 ~ R
		# y: 0 ~ R, l_center=(L-1)/2 ~ 0.0, L-1 ~ -R
		y_spatial = Dx*(-y_index + l_center)
		x_spatial = Dx*(x_index - l_center)
	
		return [y_spatial,x_spatial]
	# end def for find_spatial 
	
	############################
	# Declare and find the x_indexes of the dipole point charges
	l_left_Q = -99
	l_right_Q = -99
	for l_search in range(L-1):
		left_x_bound = find_spatial(l_center, l_search)[1] - 0.5*Dx
		if( left_x_bound <= -a/2 and -a/2 < left_x_bound + Dx):
			l_left_Q = l_search
		if( left_x_bound <= a/2 and a/2 < left_x_bound + Dx):
			l_right_Q = l_search
	# end search for loop
	
	# Double check that they where found successfully, exit if not
	if( l_left_Q == -99 or l_right_Q == -99 or l_left_Q == l_right_Q):
		print 'ERROR! l_left_Q or l_right_Q not found! Exiting!'
		sys.exit()
	
	########################################################
	# Define a function to make sure we don't overwrite boundary conditions
	# ie don't change for r > R or at the point charges 
	def write_allowed(y_index,x_index):
		# Check it's not a point charge
		if(y_index == l_center and (x_index == l_left_Q or x_index == l_right_Q)):
			return False
		# check that it's not on/beyond the R boundary
		spatial = find_spatial(y_index,x_index)
		r_squared = spatial[0]*spatial[0] + spatial[1]*spatial[1]
		if(r_squared >= R*R):
			return False
		# if we made it here it's just a normal point
		# and we can overwrite it with new data
		return True
	# end def for write_allowed
	
	########################################################
	# Define a function to make sure we don't try to read a nonexistent index 
	def get_value(V, y_index, x_index):
		# Check that the indexes are valid 
		if(y_index < 0 or L-1 < y_index or x_index < 0 or L-1 < x_index):
			# We can't return None because it won't add with floats
			# So just return 0.0, it shouldn't mess anything up
			# as these points will be beyond the R boundary condition
			return 0.0
		# if we made it here it's a valid index 
		# and we can return the value of V 
		return V[y_index][x_index]
	# end def for get_value 
	
	############################
	# Declare our ndarray, structured as:
	# np.zeros((L_y, L_x)), np.zeros((num_rows, num_columns))
	# V[y_index,x_index] = V[row_index, column_index]
	V = np.zeros((L, L))
	
	# Set the initial conditions, ie add point charges
	# R boundary already set to be zero and held fixed by write_allowed
	V[l_center][l_left_Q] = -Q_over_epsilon0
	V[l_center][l_right_Q] = Q_over_epsilon0
	
	# Print the initial V if debugging
	if(debugging): plot_V(' Initial', 'initial', V, 0, 10, './output/poisson_dipole', L, Dx, sim_method, halt_method, epsilon, fixed_accuracy)

	########################################################
	# Define a function to perform a sweep through the whole ndarray
	# using the Jacobi relaxation algorithm
	def jacobi_sweep(m_V):
		V_old = m_V.copy() # save the current values
	
		tol_sum = 0.0 # to find average abs tolerance, numerator of average
		sim_points = 0 # number of points that got updated, denominator of average
		# we need to define this because at the start only the neighbors of the dipole
		# points will change, and we don't want 0 -> 0 far away to count as being convergent,
		# ie 0 contribution when it really wasn't simulated at all!
	
		# loop over the array
		for l_y in range(L-1):
			for l_x in range(L-1):
				# Check if we can edit this cell
				if(write_allowed(l_y,l_x)):
					V_old_value = m_V[l_y][l_x] # save the old value for tol_sum calculation (not a necessary step, but in prep for GS method)
					# Update from V_old, Jacobi method = average nearest neighbors
					m_V[l_y][l_x] = 0.25*( get_value(V_old,l_y,l_x-1) + get_value(V_old,l_y-1,l_x) + get_value(V_old,l_y,l_x+1) + get_value(V_old,l_y+1,l_x) )
					# Compute the contribution to the tolerance for this nonboundary cell
					tol_sum += abs(V_old_value - m_V[l_y][l_x])
					sim_points += 1

		tol_sweep = tol_sum/sim_points
		return [m_V, tol_sweep]
	
	# end def for jacobi_sweep
	
	############################
	# run sweeps until convergence criteria is meet
	n_sweep = 0 # sweep number
	sweep = [V, 2*epsilon]

	if(halt_method == 'epsilon'):
		while (sweep[1] >= epsilon or n_sweep < l_center+1):
			# Force it to run for at least l_center+1 # of times so that
			# there is a chance for the origin to propagate out to the edge
			# otherwise wait till the tolerance is low enough
			if(sim_method == 'jacobi'):
				sweep = jacobi_sweep(sweep[0])
			#if(sim_method == 'SOR'):
				# TODO
			else:
				print 'ERROR! Unknown sim_method! Exiting!'
				sys.exit()
			n_sweep += 1
			if(debugging): print 'Sweep number = %5d, tol_current = %.5E' % (n_sweep, sweep[1])
	#elif(halt_method == 'fixed_accuracy'):
		# TODO
	else:
		print 'ERROR! Unknown halt_method! Exiting!'
		sys.exit()

	print 'run_sim completed'
	return [sweep[0], n_sweep] # TODO keep returning the info we need 
# end def for run_sim



########################################################
# Define a function to plot V 
def plot_V(optional_title, fname, m_V, n_sweep, n_contours, m_path, L, Dx, sim_method, halt_method, epsilon, fixed_accuracy):
	print 'Beginning plot_V:'	
	
	# Define arrays to make meshgrid 
	y = np.linspace(R,-R,L)
	x = np.linspace(-R,R,L)
	X,Y = np.meshgrid(x,y)

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title('$V\left(x,y\\right)$'+optional_title)
        ax.set_xlabel('$x$ [m]')
        ax.set_ylabel('$y$ [m]')

	# I like the default inwards ticks
	#mpl.rcParams['xtick.direction'] = 'out'
	#mpl.rcParams['ytick.direction'] = 'out'

	# Make a contour plot
	CS = ax.contour(X,Y,m_V,n_contours, antialiased=True)

	# Set the contour labels
	ax.clabel(CS, inline=1, fontsize=10, fmt='%2.2f')

	# Set the color bar
	CB = plt.colorbar(CS, shrink=0.8, extend='both', filled=True, cmap='inferno', ax=ax, label='$V\left(x,y\\right)$ [V]')

	# If it's the initial/diagnostics plot, draw vertical dashed lines at +-a
	if fname == 'initial':
		ax.axvline(x=a/2, ls = 'dashed', label='$a/2$', c='red')
		ax.axvline(x=-a/2, ls = 'dashed', label='$-a/2$', c='blue')

	# Add a dashed circle for the R boundary condition
	boundary_circle = plt.Circle((0,0), R, ls='dashed', color='grey',fill=False)
	ax.add_artist(boundary_circle)


	# Annotate
	text = '$L =$ %3.d, $\Delta x =$ %.3f [m]' % (L, Dx)
	text += '\n$R =$ %2.1f [m], $a =$ %1.1f [m]' % (R, a)
	text += '\n$Q/\epsilon_{0} =$ %1.1f [Vm]' % (Q_over_epsilon0)
	if(sim_method == 'jacobi'): text += '\nSim. Method = Jacobi'
	if(sim_method == 'SOR'): text += '\nSim. Method = SOR'
	if(halt_method == 'epsilon'): text += '\nConv. Criteria: $\epsilon <$ %.2E' % (epsilon)
 	if(halt_method == 'fixed_accuracy'): text += '\nConv. Criteria: $A <$ %.2E' % (fixed_accuracy)
	if(optional_title != ' Initial'): text += '\n$N_{\mathrm{iter}} =$ %5.d' % (n_sweep)
	plt.figtext(0.145, 0.13, text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small')

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/V_'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	print 'V plot printed'
# end def for plot_V


########################################################
########################################################
########################################################
# Finally, actually run things!
output_path = './output/poisson_dipole'

########################################################
########################################################
# Development Runs 

############################
# TODO move these parameters to a higher level function
Dx = 0.1 # Set the spacing of the grid points

# Set the number of grid points
L = int(round(2*R/Dx))
if(L % 2 == 0): L += 1 # ensure L is odd
print 'L = %d' % L

# set the tolerance
epsilon = 0.00001 # average, absolute change sweep to sweep

sim_method = 'jacobi'
halt_method = 'epsilon'
fixed_accuracy = -99.0

# TODO print the sim values, to the plot probably

# Run the sim
m_run = run_sim(L, Dx, sim_method, halt_method, epsilon, fixed_accuracy, True)

#print m_run[0]
plot_V('', 'final', m_run[0], m_run[1], 50, output_path, L, Dx, sim_method, halt_method, epsilon, fixed_accuracy)


########################################################
########################################################
# Production Runs for paper 

if(False):
	output_path += '/plots_for_paper'



########################################################
print '\n\nDone!\n'


