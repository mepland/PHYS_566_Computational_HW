#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

########################################################
# Set fixed parameters [SI units]

Q_over_epsilon0 = 1.0 # point charge Q over epsilon0 [TODO]
d = 2.2 # dipole length [TODO]
R = 5.0 # Location of Spherical Boundary Condition [TODO]

Dx = 1.0 # Set the spacing of the grid points


########################################################
# Print out starting values
print '\nBeginning poisson_dipole.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nQ over epsilon0 = %.1f' % Q_over_epsilon0

print '\nDx = %.4f' % Dx
# TODO

print '---------------------------------------------'
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

'''
########################################################
# TODO Define a function to run the simulation
def run_sim():

	return None
# end def for run_sim
'''

############################
# TODO make L a parameter, move odd check to higher level function
# Set the number of grid points
L = int(round(2*R/Dx))
# print 'L should be = %d' % L

# ensure L is odd
if(L % 2 == 0): L += 1
# print 'Made odd, L = %d' % L
print 'L = %d' % L

# Find the index of the center/origin
l_center = (L-1)/2 # +1 to get to center, -1 as we start from 0


########################################################
# Define a function to take ndarray indexes and return spatial coordinates of the BIN CENTERS 
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
	if( left_x_bound <= -d and -d < left_x_bound + Dx):
		l_left_Q = l_search
	if( left_x_bound <= d and d < left_x_bound + Dx):
		l_right_Q = l_search
# end search for loop

# Double check that they where found successfully, exit if not
if( l_left_Q == -99 or l_right_Q == -99 or l_left_Q == l_right_Q):
	print 'ERROR! l_left_Q or l_right_Q not found!'
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
# Declare our ndarrays, structured as:
# np.zeros((L_y, L_x)), np.zeros((num_rows, num_columns))
# V[y_index,x_index] = V[row_index, column_index]
V = np.zeros((L, L))

# Set the initial conditions, ie add point charges
V[l_center][l_left_Q] = -Q_over_epsilon0
V[l_center][l_right_Q] = Q_over_epsilon0


########################################################
# Define a function to perform a sweep through the whole ndarray
# using the Jacobi relaxation algorithm
def jacobi_sweep(m_V):

	V_old = m_V # save the current values

	# loop over the array
	for l_y in range(L-1):
		for l_x in range(L-1):
			# Check if we can edit this cell
			if(write_allowed(l_y,l_x)):
				# Update from V_old, Jacobi method = average nearest neighbors
				m_V[l_y][l_x] = 0.25*( get_value(V_old,l_y,l_x-1) + get_value(V_old,l_y-1,l_x) + get_value(V_old,l_y,l_x+1) + get_value(V_old,l_y+1,l_x) )
	return m_V

# end def for jacobi_sweep

# TODO Test it out!
print V
V = jacobi_sweep(V)
print V


'''
########################################################
########################################################
########################################################
# Finally, actually run things!
output_path = './output/poisson_dipole'

########################################################
########################################################
# Development Runs 


########################################################
########################################################
# Production Runs for paper 

# run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case)
# compare_runs(run1, run1_name, run2, run2_name, title, fname) 
# energy_run(run_name, omegaD, alphaD, theta0, run_end_periodsD, sim_method, lin_case)
# res_sweep(num_runs, omegaD_percent_range, alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case)
# lyapunov_ave(alphaD_num, theta0, run_end_periodsD, sim_method)

if(True):
	output_path += '/plots_for_paper'


'''
########################################################
print '\n\nDone!\n'


'''
# Display simple V for diagnostics
# Fill the center with -99 to mark it, temporarily...
V[l_center][l_center] = -99.0
# Fill the corners with 1,2,3,4 to mark them, temporarily...
V[0][0] = 1.0
V[0][L-1] = 2.0
V[L-1][L-1] = 3.0
V[L-1][0] = 4.0
# Fill the charges with +-7 to mark them, temporarily...
V[l_center][l_left_Q] = -7.0
V[l_center][l_right_Q] = 7.0
print '1.0 point = (y,x) = (%.3f,%.3f)' % (find_spatial(0,0)[0], find_spatial(0,0)[1])
print '2.0 point = (y,x) = (%.3f,%.3f)' % (find_spatial(0,L-1)[0], find_spatial(0,L-1)[1])
print '3.0 point = (y,x) = (%.3f,%.3f)' % (find_spatial(L-1,L-1)[0], find_spatial(L-1,L-1)[1])
print '4.0 point = (y,x) = (%.3f,%.3f)' % (find_spatial(L-1,0)[0], find_spatial(L-1,0)[1])
print V
'''

