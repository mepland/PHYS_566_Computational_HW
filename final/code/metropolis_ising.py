import os
import sys
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#import random
import math

########################################################
# Set fixed/global parameters

J = 1.5 # nearest neighbor interaction strength J
kB = 1.0 # Boltzmann's constant, relative to J TODO

neighborhood = 'Von Neumann'
# neighborhood = 'Moore'


########################################################
# Print out fixed values
print '\nBeginning dla.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nNN interaction strength J = %.1f ' % J
print 'Boltzmann\'s constant kB = %.1f ' % kB

print '\nNN neighborhood type = %s ' % neighborhood

print '\n---------------------------------------------'
print '---------------------------------------------\n'

########################################################
########################################################

# Set up the neighborhood of points to check
NN_list = []
if neighborhood == 'Moore':
	NN_list = [[-1,1], [0,1], [1,1], [-1,0], [1,0], [-1,-1], [0,-1], [1,-1]] # Moore neighborhood
elif neighborhood == 'Von Neumann':
	NN_list = [[0,1], [-1,0], [1,0], [0,-1]] # Von Neumann neighborhood
else:
	print 'ERROR!! Unknown neighborhood, exiting!!'
	sys.exit()


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
# Define a function to plot and fit the data
def plot(optional_title, m_path, fname, m_T, m_seed, world_grid):
	if debugging: print 'Beginning plot() for fname: '+fname	

	local_n = world_grid.shape[0] 

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title(optional_title)
        ax.set_xlabel('$x$ Index')
        ax.set_ylabel('$y$ Index')

	
        # adjust axis range
        ax.axis('scaled')
        axis_offset = 0.1*(local_n+1)
        ax.set_xlim((-axis_offset, local_n-1+axis_offset))
        ax.set_ylim((-axis_offset, local_n-1+axis_offset))

        # start list for legend entries/handles
        legend_handles = []

        Dx = 1.0 # grid spacing of the world

	firstcp1 = True
	firstcp2 = True

        # plot the world grid
	for i in range(local_n):
		for j in range(local_n):
                        if world_grid[i][j] == 1:
                                cp1 = plt.Rectangle((i-Dx/2, j-Dx/2), Dx, Dx, color='red', fill=True, label='Spin Up')
                                ax.add_artist(cp1)

                                if firstcp1:
                                        firstcp1 = False
                                        legend_handles.append(cp1)

                        if world_grid[i][j] == -1:
                                cp2 = plt.Rectangle((i-Dx/2, j-Dx/2), Dx, Dx, color='blue', fill=True, label='Spin Down')
                                ax.add_artist(cp2)

                                if firstcp2 and not firstcp1:
                                        firstcp2 = False
                                        legend_handles.append(cp2)


        # make a square on the world border
        world_border = plt.Rectangle((0-Dx/2,0-Dx/2), (local_n)*Dx, (local_n)*Dx, color='black', ls='dashed', fill=False, label='World Border')
        ax.add_artist(world_border)
        legend_handles.append(world_border)


        # draw legend
        ax.legend(handles=legend_handles, bbox_to_anchor=(1.03, 1), borderaxespad=0, loc='upper left', fontsize='x-small')

        # Annotate
        ann_text = '$n =$ %d\n$N =$ %d' % (local_n, local_n*local_n)
	ann_text += '\n$T =$ %.5g' % (m_T)
	ann_text += '\n\n$k_{\mathrm{B}} =$ %.5g' % (kB)
	ann_text += '\n$J =$ %.5g' % (J)
	ann_text += '\nRNG Seed = %d' % (m_seed)
        ann_text += '\n\nNN Neighborhood:\n'+neighborhood

        ax.text(1.0415, 0.018, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=ax.transAxes)


	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if debugging: print 'plot() completed!!!'
# end def for plot()



#######################################################
# Define a function to initialize the world grid
def initialize(n, seed):

	# Set up the numpy RNG with our seed
	np.random.seed(seed)

	# Set up the world grid with random +/-1 values
	# First randint makes 0, 1's then change 0's to -1
	world_grid = np.random.randint(2, size=(n, n))
	for i in range(n):
		for j in range(n):	
			if world_grid[i][j] == 0: world_grid[i][j] = -1

	return world_grid

# end def for initialize



#######################################################
# Define a function to perform one sweep of the world grid
def sweep(T, local_n, world_grid, NN_list = []):
	
	# local_n = world_grid.shape[0] pass as a parameter to speed up the sweeps...

	# sweep through whole grid
	for i in range(local_n):
		for j in range(local_n):

			# compute DeltaE

			NN_E_original = 0.0
			for k in range(len(NN_list)):
				NN_E_original += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%local_n ][ (j+NN_list[k][1])%local_n ]

			world_grid[i][j] = -world_grid[i][j] # test flip the spin

			# compute DeltaE for the spin flipped

			NN_E_flipped = 0.0
			for k in range(len(NN_list)):
				NN_E_flipped += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%local_n ][ (j+NN_list[k][1])%local_n ]

			DeltaE = NN_E_original - NN_E_flipped

			# if DeltaE <= 0 always keep the spin, ie do nothing as it's already been flipped
			if DeltaE > 0.0:
				# keep the spin with probability p = exp(-DeltaE/kB*T)
				p = np.exp( -DeltaE/(kB*T) )
				r = np.random.rand()

				# if p >= 1.0: print 'warning p >= 1' # Debugging

				# if r < p keep the spin, ie do nothing as it's already been flipped
				if r >= p:
					world_grid[i][j] = -world_grid[i][j] # flip the spin back to it's original position

	return world_grid

# end def for sweep


#######################################################
# Define a function to compute E of a world_grid
def E(local_n, world_grid, NN_list = []):
	
	# local_n = world_grid.shape[0]
	E = 0.0

	# sweep through whole grid
	for i in range(local_n):
		for j in range(local_n):

			# compute the contribution to E from this spin
			for k in range(len(NN_list)):
				E += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%local_n ][ (j+NN_list[k][1])%local_n ]

	return E
# end def for E


#######################################################
# Define a function to compute M of a world_grid
def M(local_n, world_grid):
	
	# local_n = world_grid.shape[0]
	M = 0.0

	# sweep through whole grid
	for i in range(local_n):
		for j in range(local_n):
				M += world_grid[i][j]

	return M/(local_n*local_n)
# end def for M

########################################################
########################################################
########################################################
# Finally, actually run things!

########################################################
########################################################
# Development Runs 

if(True):
	output_path = '../output/dev'
	debugging = True

	seed = 5
	test_n = 20
	T = 0.00000000000000000000001

	world_grid = initialize(test_n, seed)

	plot('initial', output_path, 'test_initial', T, seed, world_grid)

	halt_condition = 0.001

	num_history = 10
	history = np.linspace(99.0, 99.0, num_history)

	old_conv_var = 99.0
	new_conv_var = 99.0

	sweep_number = 0

	while np.mean(history) > halt_condition:

		world_grid = sweep(T, test_n, world_grid, NN_list)

		new_conv_var = E(test_n, world_grid, NN_list)

		history[sweep_number%num_history] = abs(old_conv_var - new_conv_var)/new_conv_var

		print '%d, DeltaE = %.5f, E = %.5f, M = %.5f, convergence = %.5f' % (sweep_number, old_conv_var - new_conv_var, new_conv_var, M(test_n, world_grid), np.mean(history) )
		old_conv_var = new_conv_var

		sweep_number += 1


	plot('converged', output_path, 'test_converged', T, seed, world_grid)



########################################################
########################################################
# Production Runs for paper 

if(False):
	top_output_path = '../output/plots_for_paper'
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

########################################################
print '\n\nDone!\n'


