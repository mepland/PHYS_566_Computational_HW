import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#import random
import math

import time
start_time = time.time()

import sweepMod

########################################################
# Set fixed/global parameters

J = 1.5 # nearest neighbor interaction strength J
kB = 1.0 # Boltzmann's constant, relative to J

nh = ''
NN_type = 0
neighborhood = 'Von Neumann'
# neighborhood = 'Moore'


########################################################
# Print out fixed values
print '\nBeginning metropolis_ising.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nNN interaction strength J = %.1f [J]' % J
print 'Boltzmann\'s constant kB = %.1f [J/K]' % kB

print '\nNN neighborhood type = %s' % neighborhood

print '\n---------------------------------------------'
print '---------------------------------------------\n'

########################################################
########################################################

# Set up the neighborhood of points to check, helper variables
if neighborhood == 'Von Neumann':
	NN_type = 0
	nh = 'von_neumann'
elif neighborhood == 'Moore':
	NN_type = 1
	nh = 'moore'
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



########################################################
# Define a function to plot C vs T
def plot_CT(optional_title, m_path, fname, C_list, T_list, m_n, sweep_upper_limit, halt_percent_change, m_seed):
	if debugging: print 'Beginning plot_CT() for fname: '+fname

	# Setup the arrays
	C_array = np.array(C_list)
	T_array = np.array(T_list)

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title(optional_title)
        ax.set_xlabel('$T$ [K]')
        ax.set_ylabel('$C$ [J/K]')

        # start list for legend entries/handles
        legend_handles = []

	# Plot the data
	points = ax.scatter(T_array, C_array, marker='o', label='C', c='blue')
	legend_handles.append(points)
	

	# adjust axis range
	x1_auto,x2_auto,y1_auto,y2_auto = ax.axis()
	ax.set_xlim(min(T_list)-0.5, max(T_list)+0.5)
	ax.set_ylim(0.0, y2_auto)

	
	# draw legend
	# ax.legend(handles=legend_handles, bbox_to_anchor=(0.98, 0.98), borderaxespad=0, loc='upper right', fontsize='x-small')
	
        # Annotate
	ann_text = '$n =$ %d, $N =$ %d' % (m_n, m_n*m_n)
	ann_text += '\nMax # Sweeps $=$ %g' % (sweep_upper_limit)	
	ann_text += '\n$\langle\Delta E / E \\rangle <$ %.3f' % (halt_percent_change)
	ann_text += '\n\n$J =$ %.5g [J]' % (J)
	ann_text += '\n$k_{\mathrm{B}} =$ %.4g [J/K]' % (kB)
	ann_text += '\nRNG Seed = %d' % (m_seed)
        ann_text += '\nNN Neighborhood:\n'+neighborhood

		
	ax.text(0.025, 0.72, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=ax.transAxes)

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if debugging: print 'plot_CT() completed!!!'
# end def for plot_CT()

########################################################
# Define a function to plot C max/N vs n
def plot_Cmax_vs_n(optional_title, m_path, fname, C_array, n_array, sweep_upper_limit, halt_percent_change, perform_fit, fit_upper_limit):
	if debugging: print '\n\n---------------------------------------------\n---------------------------------------------\nBeginning Cmax_vs_n for fname: '+fname

	# Setup the arrays

	C_over_N_fit = []
	C_over_N_nofit = []
	n_fit = []
	n_nofit = []


       	if fit_upper_limit is None:
		fit_upper_limit = 2.0*np.amax( n_array )

	for i in range(n_array.size):
		if n_array[i] < fit_upper_limit:
			C_over_N_fit.append( C_array[i]/(n_array[i]*n_array[i]) )
			n_fit.append(n_array[i])
		else:
			C_over_N_nofit.append( C_array[i]/(n_array[i]*n_array[i]) )
			n_nofit.append(n_array[i])

	C_over_N_fit_array = np.array(C_over_N_fit)
	C_over_N_nofit_array = np.array(C_over_N_nofit)
	log_n_fit_array = np.log(np.array(n_fit))
	log_n_nofit_array = np.log(np.array(n_nofit))


	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title(optional_title)
        ax.set_xlabel('$\log(n)$')
        ax.set_ylabel('$C/N$ [J/K]')


        # start list for legend entries/handles
        legend_handles = []

	# Plot the data
	fit_points = ax.scatter(log_n_fit_array, C_over_N_fit_array, marker='o', label='$C/N$', c='blue')
	legend_handles.append(fit_points)
	
	nofit_points = ax.scatter(log_n_nofit_array, C_over_N_nofit_array, marker='o', label=None, c='blue')

	# adjust axis range
	x2_scale_factor = 1.3
	y2_scale_factor = 1.0

	x1_auto,x2_auto,y1_auto,y2_auto = ax.axis()
	ax.set_xlim(0.0, x2_scale_factor*x2_auto)
	ax.set_ylim(0.0, y2_scale_factor*y2_auto)

	if(perform_fit):
	
		# Fitting 
		########################################################
	
	
		########################################################
		# Define the linear fit function, with offset
		def linear_fit_function(data, offset_fit, slope_fit):
			return offset_fit + slope_fit*data
		# end def linear_fit_function


		# actually perform the fits
		# op_par = optimal parameters, covar_matrix has covariance but no errors on plot so it's incorrect...
	
		m_p0 = [0.0, 1.0]	
		# m_p0 = [0.0, 3.6]

		fit_status = True		
		maxfev=m_maxfev = 2000	
		fit_text = ''
			
		try:
			op_par, covar_matrix = curve_fit(linear_fit_function, log_n_fit_array, C_over_N_fit_array, p0=m_p0, maxfev=m_maxfev)
		except RuntimeError:
			print sys.exc_info()[1]
			print 'curve_fit failed, continuing...'
			fit_status = False
		
		# plot the fit
		if(fit_status):
			# make nice x array for fit plot
			fit_x = np.linspace(min(0.1, 0.00001*abs(x2_auto)), x2_scale_factor*x2_auto, 1000)

			fit_line, = ax.plot(fit_x, linear_fit_function(fit_x, *op_par), ls='dashed', label='Linear Fit', c="black")
			legend_handles.append(fit_line)
	
			if fit_upper_limit != 2.0*np.amax( n_array ):
				fit_upper_limit_label = 'Fit Upper Limit %.1f' % np.log(fit_upper_limit)
				fit_upper_limit_line = ax.axvline(x=np.log(fit_upper_limit), ls = 'solid', label=fit_upper_limit_label, c='gray')
				legend_handles.append(fit_upper_limit_line)
	
		
		# Write out the fit parameters	
		fit_text = 'Fit Function: $C/N = a + b \log(n)$'
		if(fit_status):
			# fit_text += '\n$a_{\mathrm{Starting}} =$ %2.2f, $a_{\mathrm{Fit}} =$ %2.5f' % (m_p0[0], op_par[0])
			# fit_text += '\n$b_{\mathrm{Starting}} =$ %2.2f, $b_{\mathrm{Fit}} =$ %2.5f' % (m_p0[1], op_par[1])

			fit_text += '\n$a_{\mathrm{Fit}} =$ %2.5f' % (op_par[0])
			fit_text += ', $b_{\mathrm{Fit}} =$ %2.5f' % (op_par[1])

		else:
			fit_text += '\nFit Failed'
			
		# Print the fit parameters
		ax.text(0.025, 1-0.03, fit_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small', transform=ax.transAxes, va='top')
		
	
	# draw legend
	ax.legend(handles=legend_handles, bbox_to_anchor=(0.98, 0.98), borderaxespad=0, loc='upper right', fontsize='x-small')
	
        # Annotate
	ann_text = 'Max # Sweeps $=$ %g' % (sweep_upper_limit)	
	ann_text += '\n$\langle\Delta E / E \\rangle <$ %.3f' % (halt_percent_change)
	ann_text += '\n\n$J =$ %.5g [J]' % (J)
	ann_text += '\n$k_{\mathrm{B}} =$ %.4g [J/K]' % (kB)
        ann_text += '\nNN Neighborhood:\n'+neighborhood

	# Warning, make sure you can still see the n = 500 data point
	ax.text(0.025, 0.67, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small', transform=ax.transAxes)

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if debugging: print 'plot_MT() completed!!!'
# end def for plot_Cmax_vs_n()



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


'''
#######################################################
# Define a function to perform one sweep of the world grid
def sweep(T, m_n, NN_type, world_grid):

	# Set up the neighborhood of points to check
	NN_list = []
	if NN_type == 1:
		NN_list = [[-1,1], [0,1], [1,1], [-1,0], [1,0], [-1,-1], [0,-1], [1,-1]] # Moore neighborhood
	elif NN_type == 0:
		NN_list = [[1,0], [-1,0], [0,1], [0,-1]] # Von Neumann neighborhood
	else:
		print 'ERROR!! Unknown neighborhood, exiting!!'
		sys.exit()
	
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
'''

'''
#######################################################
# Define a function to compute E of a world_grid
def E(m_n, world_grid):

	# Set up the neighborhood of points to check
	NN_list = []
	if NN_type == 1:
		NN_list = [[-1,1], [0,1], [1,1], [-1,0], [1,0], [-1,-1], [0,-1], [1,-1]] # Moore neighborhood
	elif NN_type == 0:
		NN_list = [[1,0], [-1,0], [0,1], [0,-1]] # Von Neumann neighborhood
	else:
		print 'ERROR!! Unknown neighborhood, exiting!!'
		sys.exit()
		
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
'''

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
def loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, initial_seed, world_grid):

	num_history = 10
	history = np.linspace(99.0*m_n*m_n, 99.0*m_n*m_n, num_history)

	old_conv_var = 99.0*m_n*m_n
	new_conv_var = 99.0*m_n*m_n

	sweep_number = 0

	# sweep until the mean of the last num_history percent changes in the convergence variable is < halt_percent_change
	while np.mean(history) > halt_percent_change and sweep_number < sweep_upper_limit:

		# world_grid = sweep(T, m_n, NN_type, world_grid)
		sweepMod.run_sweep(T, kB, J, NN_type, initial_seed+sweep_number, world_grid)

		# new_conv_var = E(m_n, world_grid) # Use E as the convergence variable
		new_conv_var = sweepMod.find_E(J, NN_type, world_grid)

		if new_conv_var != 0.0:
			history[sweep_number%num_history] = abs((new_conv_var - old_conv_var)/new_conv_var) # store the percent change
		else:
			# if debugging2: print 'avoiding divide by zero error'
			history[sweep_number%num_history] = 10*halt_percent_change # give it a large non zero number just to kick it back up and keep it moving


#		if info: print 'Sweep #%d, E = %.2f, DeltaE = %.2f, M = %.2f, conv mean = %.5f' % (sweep_number, new_conv_var, new_conv_var - old_conv_var, M(m_n, world_grid), np.mean(history) ) # sweep info
		# if info and sweep_number > max(0, sweep_upper_limit-50): print 'Sweep #%d, E = %.2f, DeltaE = %.2f, M = %.2f, conv mean = %.5f' % (sweep_number, new_conv_var, new_conv_var - old_conv_var, M(m_n, world_grid), np.mean(history) ) # sweep info

		old_conv_var = new_conv_var

		sweep_number += 1

	if sweep_number >= sweep_upper_limit: print 'Warning sweep upper limit hit!\nOn T = %.2f' % T # TODO

	return [world_grid, sweep_number]

# end def for loop_till_conv

#######################################################
# Define a function to cool down through temps, saving M vs T for graphing
def cool_down(m_path, halt_percent_change, sweep_upper_limit, m_n, seed, temps_to_test = []):

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


		# loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, initial_seed, world_grid)	
		world_grid, sweep_number = loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, seed+(i*(sweep_upper_limit+2)), world_grid)

		plot_grid('', m_path, 'converged_'+temp, T, sweep_number, halt_percent_change, seed, world_grid)

		Ms.append( M(m_n, world_grid) )
		Ts.append( T )
		Sweeps.append( sweep_number )

		print 'T = %.2g Completed!' % T

	return [Ms, Ts, Sweeps]

# end def for cool_down

#######################################################
# Define a function to compute C for a given n, T
def C(T, num_microstates, microstate_sweep_separation, halt_percent_change, sweep_upper_limit, m_n, m_path, initial_seed, initial_world_grid):
	if info: print 'Beginning C for T = %.3f, n = %d' % (T, m_n)

	E_values = []

	# initialize and run to the first convergence	
	if info2: print '\nStarting microstate # 0'

	# loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, initial_seed, initial_world_grid)
	equalized_world_grid, sweep_number = loop_till_conv(T, halt_percent_change, sweep_upper_limit, m_n, initial_seed, initial_world_grid)
	world_grid = np.copy(equalized_world_grid)

	# E_values.append( E(m_n, world_grid) )
	E_values.append( sweepMod.find_E(J, NN_type, equalized_world_grid) )

	if debugging_plots:
		# save debugging plots
		name = 'C_T%.2f_n%d_equalized_world_grid' % (T, m_n)
		plot_grid('Equalized', m_path, name, T, None, None, -9, equalized_world_grid)

	# Run more sweeps to generate more E values
	for i in range(num_microstates-1):
		if info2: print '\nStarting microstate # %d' % (i+1)

		for j in range(microstate_sweep_separation):
			# world_grid = sweep(T, m_n, NN_type, world_grid)	
			sweepMod.run_sweep(T, kB, J, NN_type, initial_seed+1+i, world_grid)

		# E_values.append( E(m_n, world_grid) )
		E_values.append( sweepMod.find_E(J, NN_type, world_grid) )

	if debugging_plots:
		# save debugging plots
		name = 'C_T%.2f_n%d_last_sweep_world_grid' % (T, m_n)
		plot_grid('Last Sweep', m_path, name, T, None, None, -9, world_grid)


	# Compute the variance, then C
	variance_E = np.mean( np.square(E_values) ) - np.mean(E_values)**2

	C = variance_E / (kB*T*T)

	if info: print 'C for T = %.3f, n = %d Completed! C = %.4f' % (T, m_n, C)

	return [C, equalized_world_grid]
# end def for C

#######################################################
# Define a function to compute C for a given n, T
def CT_for_n(m_seed, num_microstates, microstate_sweep_separation, halt_percent_change, sweep_upper_limit, m_n, m_path, temps_to_test = []):
	if debugging: print '\n---------------------------------------------\nBeginning CT_for_n with n = %d' % (m_n)

	C_values = []
	T_values = []

	initial_world_grid = initialize(m_n, m_seed)

	name = 'CT_for_n%d' % m_n

	# make sure we are going through temps high to low
	temps_to_test = sorted(temps_to_test)

	for i in range(len(temps_to_test)):

		T = temps_to_test[len(temps_to_test)-1-i]

	# 	print '\nStarting T = %.2g' % T

		if i == 0:
			starting_world_grid = initial_world_grid
		else:
			starting_world_grid = equalized_world_grid

		# C(T, num_microstates, microstate_sweep_separation, halt_percent_change, sweep_upper_limit, m_n, m_path, initial_seed, initial_world_grid)
		tmp_C, equalized_world_grid = C(T, num_microstates, microstate_sweep_separation, halt_percent_change, sweep_upper_limit, m_n, m_path+'/'+name+'_debug', m_seed+1+i, starting_world_grid)
	

		C_values.append( tmp_C )
		T_values.append( T )

	plot_CT('', m_path, name, C_values, T_values, m_n, sweep_upper_limit, halt_percent_change, m_seed)

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

if(False):
	output_path = '../output/dev_'+nh
	debugging = True
	debugging2 = False
	debugging_plots = False
	info = False
	info2 = False


	# C module testing
	T = 1.0
	seed = 7
	n = 500

	initial_world_grid = initialize(n, seed)
	world_grid = np.copy(initial_world_grid)

	for i in range(15):
		print 'On Sweep #%d' % i
		seed += 1
		sweepMod.run_sweep(T, kB, J, NN_type, seed, world_grid)

	t1 = time.time()
	print 'Python E = %.5f' % ( E(n, world_grid) )
	print('Run Time: %s seconds' % (time.time() - t1)) 

	t2 = time.time()
	print 'Compiled C E = %.5f' % ( sweepMod.find_E(J, NN_type, world_grid) )  
	print('Run Time: %s seconds' % (time.time() - t2)) 

	# plot_grid('Initial', output_path, 'initial', T, None, None, -9, initial_world_grid)
	# plot_grid('Final', output_path, 'final', T, None, None, seed, world_grid)


########################################################
########################################################
# Production Runs for paper 

if(True):
	top_output_path = '../output/plots_for_paper_'+nh
	debugging = True
	debugging_plots = False # Good to see, but slows the process down ridiculously 
	debugging2 = True
	info = False
	info2 = False

        # Part a
        ########################################################
	if(True):
	        print '\nPart a:'
	        output_path = top_output_path+'/part_a'
	
		sweep_upper_limit = 20000
		n = 50
		seed = 7
		halt_percent_change = 0.01
	
		temps_to_test = [0.1, 0.5, 1, 1.5, 2, 2.25, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.7, 3.8, 3.9, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9]
	
		generate_fresh = False
	
		if os.path.isfile(output_path+'/M_array.npy') and os.path.isfile(output_path+'/T_array.npy') and os.path.isfile(output_path+'/Sweeps_array.npy') and not generate_fresh:
			print 'Loading saved arrays'
		else:
			print 'Generating M, T and Sweeps arrays'
	
			MTSweeps = cool_down(output_path, halt_percent_change, sweep_upper_limit, n, seed, temps_to_test)
	
			M_array = np.array(MTSweeps[0])
			np.save(output_path+'/M_array', M_array)
	
			T_array = np.array(MTSweeps[1])
			np.save(output_path+'/T_array', T_array)
	
			Sweeps_array = np.array(MTSweeps[2])
			np.save(output_path+'/Sweeps_array', Sweeps_array)
	
	
		M_array = np.load(output_path+'/M_array.npy')
		T_array = np.load(output_path+'/T_array.npy')
		Sweeps_array = np.load(output_path+'/Sweeps_array.npy')
	
		plot_MT('', output_path, 'M_vs_T', M_array, T_array, Sweeps_array, n, sweep_upper_limit, halt_percent_change, seed, 3.5)

	
	# Part b
        ########################################################
	if(False):
	        print '\nPart b:'
		output_path = top_output_path+'/part_b'


		sweep_upper_limit = 10000
	
		seed = 7
		halt_percent_change = 0.01
	
		num_microstates = 100
		microstate_sweep_separation = 10
	
		temps_to_test = [2.0, 2.5, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 5.0]	
		n_to_test = [5, 10, 20, 30, 40, 50, 75, 100, 200, 500]
	
		generate_fresh = False
	
		C_max_list = []
	
		if os.path.isfile(output_path+'/C_array.npy') and os.path.isfile(output_path+'/n_array.npy') and not generate_fresh:
			print 'Loading saved arrays'
		else:
			print 'Generating C and n arrays'
	
			for n in n_to_test:
	
				Ts, Cs = CT_for_n(seed, num_microstates, microstate_sweep_separation, halt_percent_change, sweep_upper_limit, n, output_path, temps_to_test)
	
				C_max_list.append(max(Cs))
	
				print '\nC max = %.5f' % (max(Cs))
	
			C_array = np.array(C_max_list)
			np.save(output_path+'/C_array', C_array)
	
			n_array = np.array(n_to_test)
			np.save(output_path+'/n_array', n_array)
	
	
		C_array = np.load(output_path+'/C_array.npy')
		n_array = np.load(output_path+'/n_array.npy')
	
		plot_Cmax_vs_n('', output_path, 'Cmax_over_N_vs_n', C_array, n_array, sweep_upper_limit, halt_percent_change, True, 60.0)


########################################################
print '\n\nDone!\n'

print('Run Time: %s seconds' % (time.time() - start_time)) 
