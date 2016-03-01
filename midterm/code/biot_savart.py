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

########################################################
# Set fixed/global parameters [SI units]

L = 1.0 # side length [m]
mu0I = 8.0 # mu0*I [1257 T nm]
m_T_con_factor = 1257.0*10**-9 # TODO use this for axes scaling or not?

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
	if( cross_product[0] == 0 and cross_product[1] == 0 and cross_product[2] == 0):
		if(debugging): print 'cross_product = 0'

		return np.array([0,0,0])
	else:
		return (mu0I/(4.0*np.pi*s_mag**3))*cross_product
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
########################################################
########################################################
# Finally, actually run things!

########################################################
########################################################
# Development Runs 
output_path = './output/development/biot_savart'

debugging = True

n = 100
Dx = L/float(n)
r = np.array([0.0,0.0,0.1])

print loop_integral(r, Dx, n)

# TODO figure out the units for mu0 I / B

# write plotting function, run it for the desired points/lines


########################################################
########################################################
# Production Runs for paper 

if(False):
	top_output_path = './output/plots_for_paper/biot_savart'
	debugging = False

        # Part a
        ########################################################
        print '\nPart a:'
        output_path = top_output_path+'/part_a'

	# TODO

        # Part b
        ########################################################
        print '\nPart b:'
        output_path = top_output_path+'/part_b'

	# TODO 

        # Part c
        ########################################################
        print '\nPart c:'
        output_path = top_output_path+'/part_c'

	# TODO

########################################################
print '\n\nDone!\n'

'''
########################################################
# Define a function to plot
def plot_V(optional_title, fname, n_contours, m_path, run=[]):
	if(debugging): print 'Beginning plot_V for fname: '+fname	

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title('$V\left(x,y\\right)$'+optional_title)
        ax.set_xlabel('$x$ [m]')
        ax.set_ylabel('$y$ [m]')


	# Annotate
	ann_text = '$L =$ %3.d, $\Delta x =$ %.3f [m]' % (L, Dx)
	plt.figtext(0.145, 0.13, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small')

	# set axis
	ax_max = 1.1*max_spatial # Zoom out 10$
	ax.set_xlim((-ax_max,ax_max))
	ax.set_ylim((-ax_max,ax_max))

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/V_'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if(debugging): print 'V plot printed'
# end def for plot_V
'''

