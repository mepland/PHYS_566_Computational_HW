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

########################################################
# Print out starting values
print '\nBeginning chaotic_pendulum.py simulation'
print '\nFixed Parameters are:'
print '---------------------------------------------'

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


########################################################
# Define a function to run the simulation
def run_sim():

	return None
# end def for run_sim



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



########################################################
print '\n\nDone!\n'


