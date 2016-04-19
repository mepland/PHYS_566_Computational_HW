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
kB = 1.0 # Boltzmann's constant, relative to J

nh = ''
neighborhood = 'Von Neumann'
# neighborhood = 'Moore'


########################################################
# Print out fixed values
print '\nBeginning dla.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nNN interaction strength J = %.1f [J]' % J
print 'Boltzmann\'s constant kB = %.1f [J/K]' % kB

print '\nNN neighborhood type = %s' % neighborhood

print '\n---------------------------------------------'
print '---------------------------------------------\n'

########################################################
########################################################

# Set up the neighborhood of points to check
NN_list = []
if neighborhood == 'Moore':
	NN_list = [[-1,1], [0,1], [1,1], [-1,0], [1,0], [-1,-1], [0,-1], [1,-1]] # Moore neighborhood
	nh = 'moore'
elif neighborhood == 'Von Neumann':
	NN_list = [[1,0], [-1,0], [0,1], [0,-1]] # Von Neumann neighborhood
	nh = 'von_neumann'
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
# Define a function to plot the world grid
def plot_grid(optional_title, m_path, fname, m_T, m_sweep_number, halt_percent_change, m_seed, world_grid):
	if debugging: print 'Beginning plot_grid() for fname: '+fname	

	m_n = world_grid.shape[0] 

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title(optional_title)
        ax.set_xlabel('$x$ Index')
        ax.set_ylabel('$y$ Index')

	
        # adjust axis range
        ax.axis('scaled')
        axis_offset = 0.1*(m_n+1)
        ax.set_xlim((-axis_offset, m_n-1+axis_offset))
        ax.set_ylim((-axis_offset, m_n-1+axis_offset))

        # start list for legend entries/handles
        legend_handles = []

        Dx = 1.0 # grid spacing of the world

	firstcp1 = True
	firstcp2 = True

        # plot the world grid
	for i in range(m_n):
		for j in range(m_n):
                        if world_grid[i][j] == 1:
                                cp1 = plt.Rectangle((i-Dx/2, j-Dx/2), Dx, Dx, color='orange', fill=True, label='Spin Up')
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
        world_border = plt.Rectangle((0-Dx/2,0-Dx/2), (m_n)*Dx, (m_n)*Dx, color='black', ls='dashed', fill=False, label='World Border')
        ax.add_artist(world_border)
        legend_handles.append(world_border)


        # draw legend
        ax.legend(handles=legend_handles, bbox_to_anchor=(1.025, 1), borderaxespad=0, loc='upper left', fontsize='x-small')

        # Annotate
       	ann_text = '$T =$ %.5g [K]' % (m_T)
       	if not m_sweep_number is None:
		ann_text += '\n# Sweeps $=$ %g' % (m_sweep_number)
		if m_sweep_number >= sweep_upper_limit: ann_text += '\nHit Limit'
		ann_text += '\n$\langle\Delta E / E \\rangle <$ %.3f' % (halt_percent_change)

	ann_text += '\n\n$n =$ %d\n$N =$ %d' % (m_n, m_n*m_n)

	ann_text += '\n\n$J =$ %.5g [J]' % (J)
	ann_text += '\n$k_{\mathrm{B}} =$ %.4g [J/K]' % (kB)
	ann_text += '\nRNG Seed = %d' % (m_seed)
        ann_text += '\nNN Neighborhood:\n'+neighborhood

        ax.text(1.04, 0.018, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=ax.transAxes)


	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if debugging: print 'plot_grid() completed!!!'
# end def for plot_grid()



########################################################
# Define a function to plot M vs T
def plot_MT(optional_title, m_path, fname, M_array, T_array, Sweeps_array, m_n, sweep_upper_limit, halt_percent_change, m_seed, Tc):
	if debugging: print 'Beginning plot_MT() for fname: '+fname

	# Setup the arrays

	M_conv = []
	M_timedout = []
	T_conv = []
	T_timedout = []

	for i in range(Sweeps_array.size):
		if Sweeps_array[i] < sweep_upper_limit:
			M_conv.append(M_array[i])
			T_conv.append(T_array[i])
		else:
			M_timedout.append(M_array[i])
			T_timedout.append(T_array[i])

	M_conv_array = np.array(M_conv)
	M_timedout_array = np.array(M_timedout)
	T_conv_array = np.array(T_conv)
	T_timedout_array = np.array(T_timedout)


	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title(optional_title)
        ax.set_xlabel('$T$ [K]')
        ax.set_ylabel('$M = N\langle s\\rangle$ [$s_{i}$]')

        # start list for legend entries/handles
        legend_handles = []

	# Plot the data
	conv_points = ax.scatter(T_conv, M_conv, marker='o', label='Converged', c='blue')
	legend_handles.append(conv_points)
	
	timedout_points = ax.scatter(T_timedout, M_timedout, marker='s', label='Timed Out', c='Green')
	legend_handles.append(timedout_points)
	
	# Plot the Tc line, if given
       	if not Tc is None:
		Tc_label = '$T_{C}$ = %.2f [K]' % Tc
		Tc_line = ax.axvline(x=Tc, ls = 'dashed', label=Tc_label, c='gray')
		legend_handles.append(Tc_line)

	# adjust axis range
	x1_auto,x2_auto,y1_auto,y2_auto = ax.axis()
	ax.set_xlim(0.0, x2_auto)
	ax.set_ylim(-max(abs(y1_auto), abs(y2_auto)), max(abs(y1_auto), abs(y2_auto)) )

	
	# draw legend
	ax.legend(handles=legend_handles, bbox_to_anchor=(0.98, 0.98), borderaxespad=0, loc='upper right', fontsize='x-small')
	
        # Annotate
	ann_text = '$n =$ %d, $N =$ %d' % (m_n, m_n*m_n)
	ann_text += '\nMax # Sweeps $=$ %g' % (sweep_upper_limit)	
	ann_text += '\n$\langle\Delta E / E \\rangle <$ %.3f' % (halt_percent_change)
	ann_text += '\n\n$J =$ %.5g [J]' % (J)
	ann_text += '\n$k_{\mathrm{B}} =$ %.4g [J/K]' % (kB)
	ann_text += '\nRNG Seed = %d' % (m_seed)
        ann_text += '\nNN Neighborhood:\n'+neighborhood

	if abs(y2_auto) < abs(y1_auto):
		# M goes positive, put ann text in upper left
	        ax.text(0.022, 0.725, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=ax.transAxes)
	else:
		# M goes negative, put ann text in lower left
	        ax.text(0.022, 0.035, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=ax.transAxes)


	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if debugging: print 'plot_MT() completed!!!'
# end def for plot_MT()



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
def sweep(T, m_n, world_grid, NN_list = []):
	
	# m_n = world_grid.shape[0] pass as a parameter to speed up the sweeps...

	# sweep through whole grid
	for i in range(m_n):
		for j in range(m_n):

			# compute DeltaE

			E_NN_original = 0.0
			for k in range(len(NN_list)):
				E_NN_original += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%m_n ][ (j+NN_list[k][1])%m_n ]

			world_grid[i][j] = -world_grid[i][j] # test flip the spin

			# compute DeltaE for the spin flipped

			E_NN_flipped = 0.0
			for k in range(len(NN_list)):
				E_NN_flipped += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%m_n ][ (j+NN_list[k][1])%m_n ]

			DeltaE = E_NN_flipped - E_NN_original

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
def E(m_n, world_grid, NN_list = []):
	
	# m_n = world_grid.shape[0]
	E = 0.0

	# sweep through whole grid
	for i in range(m_n):
		for j in range(m_n):

			# compute the contribution to E from this spin
			for k in range(len(NN_list)):
				E += -J*world_grid[i][j]*world_grid[ (i+NN_list[k][0])%m_n ][ (j+NN_list[k][1])%m_n ]	


	return E
# end def for E


#######################################################
# Define a function to compute M of a world_grid
def M(m_n, world_grid):
	
	# m_n = world_grid.shape[0]
	M = 0.0

	# sweep through whole grid
	for i in range(m_n):
		for j in range(m_n):
				M += world_grid[i][j]

	return M/(m_n*m_n)
# end def for M

#######################################################
# Define a function to run sweeps on the world grid till convergence is met
def loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, world_grid, NN_list):

	num_history = 10
	history = np.linspace(99.0*m_n*m_n, 99.0*m_n*m_n, num_history)

	old_conv_var = 99.0*m_n*m_n
	new_conv_var = 99.0*m_n*m_n

	sweep_number = 0

	# sweep until the mean of the last num_history percent changes in the convergence variable is < halt_percent_change
	while np.mean(history) > halt_percent_change and sweep_number < sweep_upper_limit:

		world_grid = sweep(T, m_n, world_grid, NN_list)

		new_conv_var = E(m_n, world_grid, NN_list) # Use E as the convergence variable

		if new_conv_var != 0.0:
			history[sweep_number%num_history] = abs((new_conv_var - old_conv_var)/new_conv_var) # store the percent change
		else:
			# if debugging2: print 'avoiding divide by zero error'
			history[sweep_number%num_history] = 10*halt_percent_change # give it a large non zero number just to kick it back up and keep it moving


#		if info: print 'Sweep #%d, E = %.2f, DeltaE = %.2f, M = %.2f, conv mean = %.5f' % (sweep_number, new_conv_var, new_conv_var - old_conv_var, M(m_n, world_grid), np.mean(history) ) # sweep info
		if info and sweep_number > max(0, sweep_upper_limit-50): print 'Sweep #%d, E = %.2f, DeltaE = %.2f, M = %.2f, conv mean = %.5f' % (sweep_number, new_conv_var, new_conv_var - old_conv_var, M(m_n, world_grid), np.mean(history) ) # sweep info

		old_conv_var = new_conv_var

		sweep_number += 1

	if sweep_number >= sweep_upper_limit: print 'Warning sweep upper limit hit!'

	return [world_grid, sweep_number]

# end def for loop_till_conv

#######################################################
# Define a function to cool down through temps, saving M vs T for graphing
def cool_down(m_path, halt_percent_change, sweep_upper_limit, m_n, seed, temps_to_test = [], NN_list = []):

	Ms = []
	Ts = []
	Sweeps = []

	# make sure we are going through temps high to low
	temps_to_test = sorted(temps_to_test)

	for i in range(len(temps_to_test)):

		T = temps_to_test[len(temps_to_test)-1-i]

		print '\nStarting T = %.2g' % T
		temp = 'T%.2f' % T
	
		if i == 0:
			world_grid = initialize(m_n, seed)
			plot_grid('Initial', m_path, 'initial', T, None, None, seed, world_grid)

	
		world_grid, sweep_number = loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, world_grid, NN_list)

		title = ''

		if sweep_number >= sweep_upper_limit and debugging:
			title = 'Timed Out'

		plot_grid(title, m_path, 'converged_'+temp, T, sweep_number, halt_percent_change, seed, world_grid)

		Ms.append( M(m_n, world_grid) )
		Ts.append( T )
		Sweeps.append( sweep_number )

		print 'T = %.2g Completed!' % T

	return [Ms, Ts, Sweeps]

# end def for cool_down

#######################################################
# Define a function to compute C for a given n, T
def C(T, num_microstates, microstate_sweep_seperation, halt_percent_change, sweep_upper_limit, m_n, m_path, initial_world_grid, NN_list = []):
	if info: print 'Beginning C for T = %.3f, n = %d' % (T, m_n)

	E_values = []

	# initialize and run to the first convergence	
	if info: print '\nStarting microstate # 0'

	world_grid = np.copy(initial_world_grid)

	world_grid, sweep_number = loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, world_grid, NN_list)

	equalized_world_grid = np.copy(world_grid)

	E_values.append( E(m_n, world_grid, NN_list) )

	if debugging_plots:
		# save debugging plots
		name = 'C_equalized_world_grid_T%.2f_n%d' % (T, m_n)
		plot_grid('Equalized', m_path+'/debugging', name, T, None, None, -9, world_grid)

	# Run more sweeps to generate more E values
	for i in range(num_microstates-1):
		if info: print '\nStarting microstate # %d' % (i+1)

		for j in range(microstate_sweep_seperation):
			world_grid = sweep(T, m_n, world_grid, NN_list)	

		E_values.append( E(m_n, world_grid, NN_list) )

	if debugging_plots:
		# save debugging plots
		name = 'C_last_sweep_world_grid_T%.2f_n%d' % (T, m_n)
		plot_grid('Last Sweep', m_path+'/debugging', name, T, None, None, -9, world_grid)


	# Compute the variance, then C
	variance_E = np.mean( np.square(E_values) ) - np.mean(E_values)**2

	C = variance_E / (kB*T*T)

	if info: print 'C for T = %.3f, n = %d Completed! C = %.4f' % (T, m_n, C)

	return [C, equalized_world_grid]
# end def for C

#######################################################
# Define a function to compute C for a given n, T
def CT_for_n(m_seed, num_microstates, microstate_sweep_seperation, halt_percent_change, sweep_upper_limit, m_n, m_path, NN_list = [], temps_to_test = []):
	if debugging: print 'Beginning CT_for_n with n = %d' % (m_n)

	C_values = []
	T_values = []

	initial_world_grid = initialize(m_n, m_seed)

	# make sure we are going through temps high to low
	temps_to_test = sorted(temps_to_test)

	for i in range(len(temps_to_test)):

		T = temps_to_test[len(temps_to_test)-1-i]

		print '\nStarting T = %.2g' % T
#		temp = 'T%.2f' % T

		if i != 0:
			starting_world_grid = equalized_world_grid
		else:
			starting_world_grid = initial_world_grid

		tmp_C, equalized_world_grid = C(T, num_microstates, microstate_sweep_seperation, halt_percent_change, sweep_upper_limit, m_n, m_path, starting_world_grid, NN_list)

		C_values.append( tmp_C )
		T_values.append( T )


	if debugging: print 'CT_for_n with n = %d Completed!' % (m_n)

	return [T_values, C_values]

# end def for CT_for_n

########################################################
########################################################
########################################################
# Finally, actually run things!

########################################################
########################################################
# Development Runs 

if(True):
	output_path = '../output/dev_'+nh
	debugging = True
	debugging2 = True
	debugging_plots = True
	info = True

	sweep_upper_limit = 2000
	m_n = 20
	m_seed = 7
	halt_percent_change = 0.01

	T = 2.0
	num_microstates = 100
	microstate_sweep_seperation = 10

	m_path = output_path

	temps_to_test = [1.0, 4.0, 8.0]

	
	Ts, Cs = CT_for_n(m_seed, num_microstates, microstate_sweep_seperation, halt_percent_change, sweep_upper_limit, m_n, m_path, NN_list, temps_to_test)



########################################################
########################################################
# Production Runs for paper 

if(False):
	top_output_path = '../output/plots_for_paper_'+nh
	debugging = True
	debugging_plots = False
	debugging2 = False
	info = True

        # Part a
        ########################################################
        print '\nPart a:'
        output_path = top_output_path+'/part_a'

	sweep_upper_limit = 5000
	n = 50
	seed = 7
	halt_percent_change = 0.01

	temps_to_test = [0.1, 0.5, 1, 1.5, 2, 2.25, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.7, 3.8, 3.9, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9]

	generate_fresh = True

	if os.path.isfile(output_path+'/M_array.npy') and os.path.isfile(output_path+'/T_array.npy') and os.path.isfile(output_path+'/Sweeps_array.npy') and not generate_fresh:
		print 'Loading saved arrays'
	else:
		print 'Generating M, T and Sweeps arrays'

		MTSweeps = cool_down(output_path, halt_percent_change, sweep_upper_limit, n, seed, temps_to_test, NN_list)

		M_array = np.array(MTSweeps[0])
		np.save(output_path+'/M_array', M_array)

		T_array = np.array(MTSweeps[1])
		np.save(output_path+'/T_array', T_array)

		Sweeps_array = np.array(MTSweeps[2])
		np.save(output_path+'/Sweeps_array', Sweeps_array)


	M_array = np.load(output_path+'/M_array.npy')
	T_array = np.load(output_path+'/T_array.npy')
	Sweeps_array = np.load(output_path+'/Sweeps_array.npy')

	plot_MT('', output_path, 'M_vs_T', M_array, T_array, Sweeps_array, n, sweep_upper_limit, halt_percent_change, seed, 3.46)


	# Part b
        ########################################################
        print '\nPart b:'
        output_path = top_output_path+'/part_b'

	# TODO

########################################################
print '\n\nDone!\n'


