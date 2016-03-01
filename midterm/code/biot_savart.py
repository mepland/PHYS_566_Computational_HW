#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# from numpy import linalg as LA
# Take the norm yourself, don't use LA.norm as it's slower for simple 1D vectors
import math
from matplotlib.ticker import FormatStrFormatter

########################################################
# Set fixed/global parameters [SI units]

L = 1.0 # side length [m]
mu0I = 8.0 # mu0*I [1257 T nm]

# Unit conversion
m_T_conv_factor = 1257.0*10**-9
B_units = 'T'

########################################################
# Print out fixed values
print '\nBeginning biot_savart.py.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nSide Length L = %.1f [m]' % L
print '"Current" mu0*I = %.1f [1257 T nm]' % mu0I

print '\n---------------------------------------------'
print '---------------------------------------------\n'

########################################################
########################################################

# define some useful vectors, I wish I could make them constant....

x = np.array([1.0,0.0,0.0])
y = np.array([0.0,1.0,0.0])
z = np.array([0.0,0.0,1.0])

corner_A = np.array([ L/2.0, -L/2.0, 0.0])
corner_B = np.array([ L/2.0,  L/2.0, 0.0])
corner_C = np.array([-L/2.0,  L/2.0, 0.0])
corner_D = np.array([-L/2.0, -L/2.0, 0.0])

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
# Define the integrand function
def dB(l_hat, r, r_prime):

	# debugging flagged safety check that l_hat is a unit vector
	if(debugging):
		if( math.sqrt(l_hat[0]**2 + l_hat[1]**2 + l_hat[2]**2) != 1.0):
			print 'ERROR! l_hat is not of unit length! Exiting!'
			sys.exit()

	s = r - r_prime
	s_mag = math.sqrt(s[0]**2 + s[1]**2 + s[2]**2)

	cross_product = np.cross(l_hat, s)

	# if s and our line integral curve are collinear l_hat cross s = 0
	# and we don't want to have any contribution from the integral
	# so, return 0 vector right away. This will also avoid s_mag = 0 infinity problems
	if( debugging2 and cross_product[0] == 0 and cross_product[1] == 0 and cross_product[2] == 0):
		if(debugging): print 'cross_product = 0'

		return np.array([0,0,0])
	else:
		return ((mu0I*m_T_conv_factor)/(4.0*np.pi*s_mag**3))*cross_product
# end def for dB

########################################################
# Define a line integral function to integrate over one of the squares sides
# Using the Riemann sum method
def line_integral(l_hat, l_begin, r, Dx, n):

	# debugging flagged safety check
	if(debugging):
		if( L/n != Dx):
			print 'ERROR! n and Dx are miss matched! Exiting!'
			sys.exit()

	# Do the integral
	integral = np.array([0,0,0])
	for i in range(1,n):
		# Find the midpoint vector/location
		r_prime = l_begin + l_hat*Dx*(i-0.5)
		# add it up
		integral = integral + dB(l_hat, r, r_prime)*Dx

	return integral
# end def for line_integral

########################################################
# Define a loop integral function to integrate over the whole square 
def loop_integral(r, Dx, n):

	# sum up the four sides 
	l1 = line_integral( y, corner_A, r, Dx, n)
	l2 = line_integral(-x, corner_B, r, Dx, n)
	l3 = line_integral(-y, corner_C, r, Dx, n)
	l4 = line_integral( x, corner_D, r, Dx, n)

	return l1+l2+l3+l4
# end def for loop_integral

########################################################
# Define a function to plot B along x, y or z
# fixed_point is a 3 vector, free_var_num is 0,1,2 ~ the var to cycle through
def plot_B(target_Dx, fixed_point, free_var_num, var_min, var_step1, var_trans1, var_step2, var_max, fname, m_path, optional_title, plot_theory):
	if(debugging): print 'Beginning plot_B for fname: '+fname	

	# Figure out the free and fixed variables...
	if(free_var_num != 0 and free_var_num != 1 and free_var_num != 2):
		print 'ERROR! Invalid free_var_num! Exiting!'
		sys.exit()

	fixed_vars = []
	for i in [0,1,2]:
		if(i != free_var_num): fixed_vars.append(i)

	if(free_var_num == 0):
		free_var_name = 'x'
		fixed_var1_name ='y'
		fixed_var2_name ='z'
	elif(free_var_num == 1):
		free_var_name = 'y'
		fixed_var1_name ='x'
		fixed_var2_name ='z'
	elif(free_var_num == 2):
		free_var_name = 'z'
		fixed_var1_name ='x'
		fixed_var2_name ='y'


	# Set up sim parameters
	n = int(L/target_Dx)
	Dx = (L/float(n))

	# Set up arrays to hold result
	free_var_list = []
	Bx_list = []
	By_list = []
	Bz_list = []

	# Run the loop
	free_var = var_min
	while free_var <= var_max:
		if(debugging): print 'Integrating B for %s = %.4f' % (free_var_name, free_var)
		# find the current test_r 
		test_r = fixed_point
		test_r[free_var_num] = free_var

		# sim and save data 
		B_vec = loop_integral(test_r, Dx, n)

		free_var_list.append(free_var)
		Bx_list.append(B_vec[0])
		By_list.append(B_vec[1])
		Bz_list.append(B_vec[2])

		# increment along our axis
		if(abs(free_var) <= var_trans1 - var_step2): free_var += var_step1
		else: free_var += var_step2

	if(debugging): print 'Simulation Completed, Plotting'

	# Create the ndarrays to plot
	free_var_ndarray = np.array(free_var_list)
	Bx_ndarray = np.array(Bx_list)
	By_ndarray = np.array(By_list)
	Bz_ndarray = np.array(Bz_list)

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title('$B('+free_var_name+')$'+optional_title)
        ax.set_xlabel('$'+free_var_name+'$ [m]')
        ax.set_ylabel('$B$ ['+B_units+']')

	# plot Bx, By, Bz
	ax.plot(free_var_ndarray, Bz_ndarray, ls='solid', label='$B_{z}$ Simulation', c='blue', ms = 6, marker='o')
        ax.plot(free_var_ndarray, Bx_ndarray, ls='solid', label='$B_{x}$ Simulation', c='crimson', ms = 6, marker='v')
        ax.plot(free_var_ndarray, By_ndarray, ls='solid', label='$B_{y}$ Simulation', c='lime', ms = 6, marker='^')

	ann_text_y = 0.685

	# TODO plot theory
	if(plot_theory):
		if(fixed_point[fixed_vars[0]] != 0.0 and fixed_point[fixed_vars[1]] != 0.0):
			print 'ERROR! Comparing theory apples to oranges!! Exiting!'
			sys.exit()

		########################################################
		# Define a function for the expected Bz(z=z, x=y=0) of a circular loop with the same area
		# result taken from Griffiths Ex 5.6, eq (5.41)
		def Bz_theory(z):
			R = math.sqrt((L**2)/np.pi)
 			return (mu0I*m_T_conv_factor*R**2)/( 2.0*pow((R**2 + z**2), (3.0/2.0)) )
		# end def for Bz_theory

		free_var_fine = np.linspace(var_min,var_max,200)
		ax.plot(free_var_fine, Bz_theory(free_var_fine), ls='dashed', label='$B_{z}$ Theory', c='black')
		ann_text_y = 0.65


	# Annotate
	ann_text = '$n =$ %3.d, $\Delta %s =$ %.4f [m]' % (n, free_var_name, Dx)
	ann_text += '\n$%s$ =%.1f [m], $%s$ = %.1f [m]' % (fixed_var1_name, fixed_point[fixed_vars[0]], fixed_var2_name, fixed_point[fixed_vars[1]])
	ann_text += '\n%.1f < $%s$ < %.1f [m]' % (var_min, free_var_name, var_max) 
	plt.figtext(0.697, ann_text_y, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small')

	# adjust axis
	fig.subplots_adjust(left=0.17)
	majorFormatter = FormatStrFormatter('%.1E')
	ax.get_yaxis().set_major_formatter(majorFormatter)

	# draw legend
	ax.legend(fontsize='small')

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/B_'+free_var_name+'_'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if(debugging): print 'B plot printed'
# end def for plot_B



########################################################
########################################################
########################################################
# Finally, actually run things!

########################################################
########################################################
# Development Runs 
output_path = './output/development/biot_savart'

debugging = True
debugging2 = False

# n = 200
# Dx = L/float(n)

########################################################
########################################################
# Production Runs for paper 

# plot_B(target_Dx, fixed_point, free_var_num, var_min, var_step1, var_trans1, var_step2, var_max, fname, m_path, optional_title, plot_theory)

if(True):
	top_output_path = './output/plots_for_paper/biot_savart'

	debugging = False
	debugging2 = False

	n = 300
	Dx = L/float(n)

	free_var_max = 3.5

        # Part a
        ########################################################
        print '\nPart a:'
        output_path = top_output_path+'/part_a'

	plot_B(Dx, np.array([0.0,0.0,-99.0]), 2, -free_var_max, 0.1, 2.5, 0.5, free_var_max, 'Bxyz', output_path, '', True)

        # Part b
        ########################################################
        print '\nPart b:'
        output_path = top_output_path+'/part_b'

	plot_B(Dx, np.array([-99.0,0.0,1.0]), 0, -free_var_max, 0.1, 3.0, 0.5, free_var_max, 'Bxyz', output_path, '', False)

        # Part c
        ########################################################
        print '\nPart c:'
        output_path = top_output_path+'/part_c'

	plot_B(Dx, np.array([0.5,0.0,-99.0]), 2, -free_var_max, 0.05, 2.0, 0.5, free_var_max, 'Bxyz', output_path, '', False)

########################################################
print '\n\nDone!\n'


