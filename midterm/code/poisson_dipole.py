#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
#import matplotlib as mpl # We need to import all of matplotlib just to set rcParams once...
#from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from matplotlib.ticker import LogFormatterExponent 

########################################################
# Set fixed/global parameters [SI units]

Q_over_epsilon0 = 1.0 # point charge Q over epsilon0 [Vm]

a = 0.6 # dipole separation length [m]
R = 10.0 # Location of Spherical Boundary Condition [m]

########################################################
# Print out fixed values
print '\nBeginning poisson_dipole.py'
print '\nFixed Parameters are:'
print '---------------------------------------------'

print '\nQ over epsilon0 = %.1f [Vm]' % Q_over_epsilon0
print 'Dipole Separation a = %.1f [m]' % a
print 'Circular Boundary R = %.1f [m]' % R

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

#######################################################
# Define a function to take ndarray indexes and return spatial coordinates of the BIN CENTERS 
def find_spatial(y_index,x_index, L, Dx):
	l_center = (L-1)/2 # +1 to get to center, -1 as we start from 0
	# if Dx is optimally chosen, otherwise these aren't R's...
	# X: 0 ~ - R, l_center=(L-1)/2 ~ 0.0, L-1 ~ R
	# y: 0 ~ R, l_center=(L-1)/2 ~ 0.0, L-1 ~ -R
	y_spatial = Dx*(-y_index + l_center)
	x_spatial = Dx*(x_index - l_center)
	
	return [y_spatial,x_spatial]
# end def for find_spatial 

########################################################
# Define a function to run the simulation
# sim_method, halt_method are self explanatory
# epsilon is the tolerance convergence condition/limit
# fixed_accuracy is the accuracy convergence condition/limit
# Debugging prints some extra info and makes the initial plot
def run_sim(L, Dx, sim_method, halt_method, epsilon, fixed_accuracy, extra_plots, extra_plots_fname):
	if(debugging): print 'Beginning run_sim:'	

	# Check that L is odd
	if(L % 2 == 0):
		print 'ERROR! L is even! Exiting!'
		sys.exit()

	# Find the index of the center/origin
	l_center = (L-1)/2 # +1 to get to center, -1 as we start from 0

	# Find the optimal alpha value for SOR
	alpha = 2.0/(1+(np.pi/L))
	if(alpha >=2.0): print 'WARNING ALPHA >= 2.0, SOR WILL NOT CONVERGE'
	

	# Debug the coordinate system...
	if(False):
		print '\nDebugging Coordinate system and exiting:'
		print 'R = %.2f, L*Dx = .5f, Dx = %.5f' % (R, L*Dx, Dx)
		points = [[0,0],[0,L-1],[L-1,L-1],[L-1,0]]
		for i in range(len(points)):
			print '[y_index][x_index]=[%d][%d] ~ (y,x)=(%.3f,%.3f)' %(points[i][0], points[i][1], find_spatial(points[i][0], points[i][1], L, Dx)[0], find_spatial(points[i][0], points[i][1], L, Dx)[1])
		return None

	
	############################
	# Declare and find the x_indexes of the dipole point charges
	l_left_Q = -99
	l_right_Q = -99
	for l_search in range(L-1):
		left_x_bound = find_spatial(l_center, l_search, L, Dx)[1] - 0.5*Dx
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
		spatial = find_spatial(y_index,x_index, L, Dx)
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
	# multiply by 4.0 because rho should not get the 1/4 prefactor,
	# but will otherwise the way this is coded
	# Dx^2/Dx^2 cancel out
	V[l_center][l_left_Q] = -4.0*Q_over_epsilon0
	V[l_center][l_right_Q] = 4.0*Q_over_epsilon0

	# if desired, make initial V and boundary value plots 
	if(extra_plots):
		artificial_run = [V, -99, L, Dx, alpha, sim_method, halt_method, epsilon, fixed_accuracy]
		plot_V(' Initial', 'initial'+extra_plots_fname, 10, './output/plots_for_paper_final/poisson_dipole/extra_material', artificial_run)

		# Make a boundary_V filling every allowed point with 7.7
		boundary_V = np.zeros((L,L))
		for l_y in range(L-1):
			for l_x in range(L-1):
				# Check if we can edit this cell
				if(write_allowed(l_y,l_x)): boundary_V[l_y][l_x] = 7.7
		artificial_run[0] = boundary_V
		plot_V(' Boundary', 'boundary'+extra_plots_fname, 3, './output/plots_for_paper_final/poisson_dipole/extra_material', artificial_run)


	########################################################
	# Define a function to perform a sweep through the whole ndarray
	# using the Jacobi relaxation algorithm
	def jacobi_sweep(m_V):
		V_old = m_V.copy() # save the current values

		# define variables for epsilon/tolerance halting method
		tol_sum = 0.0 # to find average abs tolerance, numerator of average
		sim_points = 0 # number of points that got updated, denominator of average
		# we need to define this because at the start only the neighbors of the dipole
		# points will change, and we don't want 0 -> 0 far away to count as being convergent,
		# ie 0 contribution when it really wasn't simulated at all!
		# Also it will not count boundary value points

		# define variables for fixed accuracy halting method
		A_sweep = 0.0

		# loop over the array
		for l_y in range(L-1):
			for l_x in range(L-1):
				# Check if we can edit this cell
				if(write_allowed(l_y,l_x)):
					V_old_value = m_V[l_y][l_x] # save the old value for tol_sum calculation (not a necessary step, but in prep for GS method)
					# Update from V_old, Jacobi method = average nearest neighbors
					m_V[l_y][l_x] = 0.25*( get_value(V_old,l_y,l_x-1) + get_value(V_old,l_y-1,l_x) + get_value(V_old,l_y,l_x+1) + get_value(V_old,l_y+1,l_x) )

					abs_difference = abs(V_old_value - m_V[l_y][l_x])

					# Compute the contribution to the tolerance for this non boundary cell
					# but only if it was updated
					if(V_old_value != m_V[l_y][l_x]):
						tol_sum += abs_difference
						sim_points += 1

					# Find the accuracy, defined as largest relative change of a single grid point from one iteration to the next
					if(A_sweep < abs_difference): A_sweep = abs_difference


		tol_sweep = tol_sum/sim_points
		return [m_V, tol_sweep, A_sweep]
	
	# end def for jacobi_sweep
	
	########################################################
	# Define a function to perform a sweep through the whole ndarray
	# using the Simultaneous Over Relaxation (SOR) algorithm
	def sor_sweep(m_V):
		# define variables for epsilon/tolerance halting method
		tol_sum = 0.0 # to find average abs tolerance, numerator of average
		sim_points = 0 # number of points that got updated, denominator of average

		# define variables for fixed accuracy halting method
		A_sweep = 0.0

		# loop over the array
		# Future work, spiral 'snake' style out from origin to not bias towards better
		# values in the lower right corner compared to the rest of the simulation region 
		for l_y in range(L-1):
			for l_x in range(L-1):
				# Check if we can edit this cell
				if(write_allowed(l_y,l_x)):
					V_old_value = m_V[l_y][l_x] # save the old value

					# Run SOR method
					V_GS = 0.25*( get_value(m_V,l_y,l_x-1) + get_value(m_V,l_y-1,l_x) + get_value(m_V,l_y,l_x+1) + get_value(m_V,l_y+1,l_x) )
					DeltaV = V_GS - V_old_value
					m_V[l_y][l_x] = alpha*DeltaV + V_old_value

					abs_difference = abs(V_old_value - m_V[l_y][l_x])

					# Compute the contribution to the tolerance for this non boundary cell
					# but only if it was updated
					if(V_old_value != m_V[l_y][l_x]):
						tol_sum += abs_difference
						sim_points += 1

					# Find the accuracy, defined as largest relative change of a single grid point from one iteration to the next
					if(A_sweep < abs_difference): A_sweep = abs_difference


		tol_sweep = tol_sum/sim_points
		return [m_V, tol_sweep, A_sweep]
	
	# end def for sor_sweep


	############################
	# run sweeps until convergence criteria is meet
	n_sweep = 0 # sweep number
	sweep = [V, -99.0, -99.0] # Store V and the convergence criteria in the sweep list

	if(halt_method == 'epsilon'):
		sweep[1] = 2*epsilon
		var_to_compare = 1
		halt_value = epsilon
	elif(halt_method == 'fixed_accuracy'):
		sweep[2] = 2*fixed_accuracy 
		var_to_compare = 2
		halt_value = fixed_accuracy
	else:
		print 'ERROR! Unknown halt_method! Exiting!'
		sys.exit()

	# Do the loop
	while (sweep[var_to_compare] >= halt_value):
		if(sim_method == 'jacobi'): sweep = jacobi_sweep(sweep[0])
		elif(sim_method == 'SOR'): sweep = sor_sweep(sweep[0])
		else:
			print 'ERROR! Unknown sim_method! Exiting!'
			sys.exit()
		n_sweep += 1
		
		if(debugging or n_sweep % 500 == 0): print 'Sweep number = %5d, tol_current = %.5E, A_current = %.5E' % (n_sweep, sweep[1], sweep[2])
	# end while loop

	if(debugging): print 'run_sim completed'
	return [sweep[0], n_sweep, L, Dx, alpha, sim_method, halt_method, epsilon, fixed_accuracy]
# end def for run_sim



########################################################
# Define a function to plot V 
def plot_V(optional_title, fname, n_contours, m_path, run=[]):
	if(debugging): print 'Beginning plot_V for fname: '+fname	

	m_V = run[0]
	n_sweep = run[1]
	L = run[2]
	Dx = run[3]
	alpha = run[4]
	sim_method = run[5]
	halt_method = run[6]
	epsilon = run[7]
	fixed_accuracy = run[8]

	# Define arrays to make meshgrid 
	l_center = (L-1)/2 # +1 to get to center, -1 as we start from 0
	max_spatial = Dx*l_center
	y = np.linspace(max_spatial,-max_spatial,L)
	x = np.linspace(-max_spatial,max_spatial,L)
	X,Y = np.meshgrid(x,y)

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title('$V(x,y)$'+optional_title)
        ax.set_xlabel('$x$ [m]')
        ax.set_ylabel('$y$ [m]')

	# I like the default inwards ticks
	#mpl.rcParams['xtick.direction'] = 'out'
	#mpl.rcParams['ytick.direction'] = 'out'

	# Make a contour plot
	if(fname == 'boundary'): n_contours = [-Q_over_epsilon0, Q_over_epsilon0, 7.7] 
	try:
		CS = ax.contour(X,Y,m_V,n_contours, antialiased=True)
	except RuntimeWarning:
		print 'ax.contour threw an warning trying to plot fname = %s, it should be fine...' % fname

	# Set the contour labels, turned off for being to messy near the charges
	# ax.clabel(CS, inline=1, fontsize=10, fmt='%2.2f')

	# Set the color bar
	CB = plt.colorbar(CS, shrink=0.8, extend='both', filled=True, cmap='inferno', ax=ax, label='$V(x,y)$ [V]')

	legend_handles = []
	# If it's the initial/diagnostics plot, draw vertical dashed lines at +-a
	if fname == 'initial':
		red_line = ax.axvline(x=a/2, ls = 'dashed', label='$a/2$', c='red')
		blue_line = ax.axvline(x=-a/2, ls = 'dashed', label='$-a/2$', c='blue')
		legend_handles.append(red_line)
		legend_handles.append(blue_line)

	# Add a dashed circle for the R boundary condition
	radius_circle = plt.Circle((0,0), R, ls='dashed', color='grey', fill=False, label='$V(R) = 0$ Boundary')
	ax.add_artist(radius_circle)
	legend_handles.append(mlines.Line2D([], [], ls='dashed', color='grey', label='$V(R) = 0$ Boundary'))

	# Add dashed square for the simulation border
	boundary_square = plt.Rectangle((-max_spatial,-max_spatial), 2*max_spatial, 2*max_spatial, angle=0.0, ls='dashed', color='black', fill=False, label='Sim. Boundary')
	ax.add_artist(boundary_square)
	legend_handles.append(mlines.Line2D([], [], ls='dashed', color='black', label='Sim. Boundary'))

	if(n_sweep != -99): 
		# Add dashed square for the max distance info can propagate in Niter steps
		max_info_prop = a + Dx*n_sweep 
		info_prop_square = plt.Rectangle((-max_info_prop,-max_info_prop), 2*max_info_prop, 2*max_info_prop, angle=0.0, ls='dotted', color='maroon', fill=False, label='Prop. Boundary')
		ax.add_artist(info_prop_square)
		legend_handles.append(mlines.Line2D([], [], ls='dotted', color='maroon', label='Prop. Boundary'))


	# Annotate
	ann_text = '$L =$ %3.d, $\Delta x =$ %.3f [m]' % (L, Dx)
	ann_text += '\n$R =$ %2.1f [m], $a =$ %1.1f [m]' % (R, a)
	ann_text += '\n$Q/\epsilon_{0} =$ %1.1f [Vm]' % (Q_over_epsilon0)
	if(fname != 'initial' and fname != 'boundary'): 
		if(sim_method == 'jacobi'): ann_text += '\nSim. Method = Jacobi'
		if(sim_method == 'SOR'): ann_text += '\nSim. Method = SOR, $\\alpha =$ %.3f' % (alpha)
		if(halt_method == 'epsilon'): ann_text += '\nConv. Criteria: $\epsilon <$ %.1E [V]' % (epsilon)
 		if(halt_method == 'fixed_accuracy'): ann_text += '\nConv. Criteria: $A <$ %.1E [V]' % (fixed_accuracy)
		ann_text += '\n$N_{\mathrm{iter}} =$ %.d' % (n_sweep)
	plt.figtext(0.145, 0.13, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small')

	# set axis
	ax_max = 1.1*max_spatial # Zoom out 10%
	ax.set_xlim((-ax_max,ax_max))
	ax.set_ylim((-ax_max,ax_max))

	# draw legend
	ax.legend(handles=legend_handles, bbox_to_anchor=(1.35, -0.114), borderaxespad=0, loc='lower right', fontsize='x-small')

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/V_'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if(debugging): print 'V plot printed'
# end def for plot_V

########################################################
# Define a function to plot V(r) along the dipoles axis 
def plot_Vr_dipole_axis(optional_title, fname, m_path, run=[]):
	if(debugging): print 'Beginning plot_Vr_dipole_axis for fname: '+fname	

	m_V = run[0]
	n_sweep = run[1]
	L = run[2]
	Dx = run[3]
	alpha = run[4]
	sim_method = run[5]
	halt_method = run[6]
	epsilon = run[7]
	fixed_accuracy = run[8]

	l_center = (L-1)/2
	max_spatial = Dx*l_center

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)
        ax.set_title('$V(r)$ along Dipole Axis'+optional_title)
        ax.set_xlabel('$r$ [m]')
        ax.set_ylabel('$V$ [V]')

	########################################################
	# Define a function for the expected V(r) 
	def Vr_dipole_axis_theory(r):
 		if(abs(r) == a/2.0):
			return None
		else:
			return (Q_over_epsilon0/(4.0*np.pi))*( (1.0/abs(r-(a/2.0))) -(1.0/abs(r+(a/2.0))) )
	# end def for Vr_dipole_axis_theory

	# loop over V[l_y][l_x] save V(r) sim and calculate V(r) theory
	r = []
	Vr = []
	Vr_theory = []

	for l_x in range(L-1):
		r.append(find_spatial(l_center, l_x, L, Dx)[1])
		Vr.append(m_V[l_center][l_x])
		Vr_theory.append(Vr_dipole_axis_theory(find_spatial(l_center, l_x, L, Dx)[1]))

	# create the ndarrays to plot
	r_ndarray = np.array(r)
	Vr_ndarray = np.array(Vr)
	Vr_theory_ndarray = np.array(Vr_theory)


	# make the plots
	ax.plot(r_ndarray, Vr_ndarray, ls='solid', label='$V(r)$ Simulation', c='black')
	ax.plot(r_ndarray, Vr_theory_ndarray, ls='dashed', label='$V(r)$ Theory', c='darkgreen')

	# Draw vertical dashed lines at +-a
	ax.axvline(x=a/2, ls = 'solid', label='$a/2$', c='red')
	ax.axvline(x=-a/2, ls = 'solid', label='$-a/2$', c='blue')

	# Draw vertical lines at +-R, max_spatial, and max_info_prop
	ax.axvline(x=R, ls = 'solid', label='$V(R) = 0$ Boundary', c='black')
	ax.axvline(x=-R, ls = 'solid', label=None, c='black')
	ax.axvline(x=max_spatial, ls = 'dashed', label='Sim. Boundary', c='grey')
	ax.axvline(x=-max_spatial, ls = 'dashed', label=None, c='grey')
	max_info_prop = a + Dx*n_sweep 
	ax.axvline(x=max_info_prop, ls='dotted', label='Prop. Boundary', c='maroon')
	ax.axvline(x=-max_info_prop, ls='dotted', label=None, c='maroon')


	# Annotate
	ann_text = '$L =$ %3.d, $\Delta x =$ %.3f [m]' % (L, Dx)
	ann_text += '\n$R =$ %2.1f [m], $a =$ %1.1f [m]' % (R, a)
	ann_text += '\n$Q/\epsilon_{0} =$ %1.1f [Vm]' % (Q_over_epsilon0)
	if(sim_method == 'jacobi'): ann_text += '\nSim. Method = Jacobi'
	if(sim_method == 'SOR'): ann_text += '\nSim. Method = SOR, $\\alpha =$ %.3f' % (alpha)
	if(halt_method == 'epsilon'): ann_text += '\nConv. Criteria: $\epsilon <$ %.1E [V]' % (epsilon)
	if(halt_method == 'fixed_accuracy'): ann_text += '\nConv. Criteria: $A <$ %.1E [V]' % (fixed_accuracy)
	ann_text += '\n$N_{\mathrm{iter}} =$ %.d' % (n_sweep)
	plt.figtext(0.145, 0.13, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small')

	# set axis
	ax_max = 1.1*max_spatial # Zoom out 10%
	ax.set_xlim((-ax_max,ax_max))
	ax.set_ylim((-ax_max,ax_max))

	# draw legend
	ax.legend(loc='lower right', fontsize='x-small')

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/Vr_dipole_axis_'+fname+'.pdf')	

	fig.clf() # Clear fig for reuse

	if(debugging): print 'V(r) dipole axis plot printed'
# end def for plot_Vr_dipole_axis

########################################################
# Define a function to plot N_iter vs epsilon 
def plot_N_vs_epsilon(L, Dx, epsilon_low, epsilon_step0, tipping_point1, epsilon_step1, tipping_point2, epsilon_step2, epsilon_high, m_path):
	print 'Beginning plot_N_vs_epsilon:'	
	
	# Set fixed parameters
	sim_method = 'jacobi'
	halt_method = 'epsilon'
	fixed_accuracy = -99.0
	num_contours = 50

	# Loop over run_sim, save epsilons and n_sweeps ~ N_iter
	epsilons = []
	n_sweeps = []

	m_epsilon = epsilon_low
	while m_epsilon <= epsilon_high:
		fname = 'V_for_epsilon_%4e' % m_epsilon
		this_run = run_sim(L, Dx, sim_method, halt_method, m_epsilon, fixed_accuracy, False, '')
		plot_V(' $N_{\mathrm{iter}}(\epsilon)$ Run', fname, num_contours, m_path+'/runs', this_run)

		print 'For epsilon = %f, Niter = %d' % (m_epsilon, this_run[1])
		epsilons.append(m_epsilon)
		n_sweeps.append(this_run[1])

		if(m_epsilon < tipping_point1):
			m_epsilon += epsilon_step0
		elif(m_epsilon < tipping_point2):
			m_epsilon += epsilon_step1
		else:
			m_epsilon += epsilon_step2

	# Extract the last two parameters, as used, from the last this_run
	L = this_run[2]
	Dx = this_run[3]

	# Create the ndarrays to plot
        epsilons_ndarray = np.array(epsilons)
	n_sweeps_ndarray = np.array(n_sweeps)

	# Set up the figure and axes
        fig = plt.figure('fig')
        ax = fig.add_subplot(111)

        ax.set_title('$N_{\mathrm{iter}}(\epsilon)$')
        ax.set_xlabel('$\epsilon$ [V]')
        ax.set_ylabel('$N_{\mathrm{iter}}$')

	# scale axis
	x_scale = 0.1
	ax.set_xlim((0.0, (1+x_scale)*max(epsilons)))
	y_scale = 0.1
	ax.set_ylim(((1-y_scale)*min(n_sweeps), (1+x_scale)*max(n_sweeps)))

	# Create the plot
	ax.scatter(epsilons_ndarray, n_sweeps_ndarray, marker='o', label='$N_{\mathrm{iter}}$', c='blue')

	# Make the plot semi log, log log doesn't work, x doesn't vary by enough
	ax.set_yscale('log')

	# Draw the legend
        ax.legend(fontsize='small')

	# Annotate
	ann_text = '$L =$ %3.d, $\Delta x =$ %.3f [m]' % (L, Dx)
	ann_text += '\n$R =$ %2.1f [m], $a =$ %1.1f [m]' % (R, a)
	ann_text += '\n$Q/\epsilon_{0} =$ %1.1f [Vm]' % (Q_over_epsilon0)
	if(sim_method == 'jacobi'): ann_text += '\nSim. Method = Jacobi'
	if(sim_method == 'SOR'): ann_text += '\nSim. Method = SOR'
	if(halt_method == 'epsilon'): ann_text += '\n$\epsilon$ Conv. Criteria'
	if(halt_method == 'fixed_accuracy'): ann_text += '\n$A$ Conv. Criteria'
	plt.figtext(0.699, 0.70, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small')

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/Niter_vs_epsilon.pdf')	

	fig.clf() # Clear fig for reuse

	print 'plot_N_vs_epsilon completed!!!'
# end def for plot_N_vs_epsilon

########################################################
# Define a function to plot N_iter vs n = L*L for fixed spatial length, variable L -> variable Dx
def plot_N_vs_n(fixed_accuracy, L_low, L_step0, tipping_point1, L_step1, L_high, m_path, jacobi_n_fit_cut):
	print 'Beginning plot_N_vs_n:'	
	
	# Set fixed parameters
	halt_method = 'fixed_accuracy'
	epsilon = -99.0
	num_contours = 50

	# loop over run_sim, save n=L*L and n_sweeps ~ n_iter
	ns_fit = []
	ns_no_fit = []
	jacobi_sweeps_fit = []
	jacobi_sweeps_no_fit = []
	sor_sweeps = []

	m_L = L_low
	while m_L <= L_high:
		if(m_L % 2 == 0): m_L += 1 # make sure L is odd, no matter what L_low, L_step are given
		m_n = m_L*m_L
		fname = 'V_for_L_%d' % m_L	

		m_Dx = 2*R/(m_L-1)

		jacobi_run = run_sim(m_L, m_Dx, 'jacobi', halt_method, epsilon, fixed_accuracy, False, '')
		plot_V(' Jacobi $N_{\mathrm{iter}}(n)$ Run', 'jacobi_'+fname, num_contours, m_path+'/jacobi_runs', jacobi_run)

		sor_run = run_sim(m_L, m_Dx, 'SOR', halt_method, epsilon, fixed_accuracy, False, '')
		plot_V(' SOR $N_{\mathrm{iter}}(n)$ Run', 'sor_'+fname, num_contours, m_path+'/sor_runs', sor_run)

		print 'For L = %3.d, Dx = %.5f [m], n = %.4E, Niter = %5.d (Jacobi), %5.d (SOR)' % (m_L, m_Dx, m_n, jacobi_run[1], sor_run[1])

		if(m_n < jacobi_n_fit_cut or jacobi_n_fit_cut == -99):
			ns_fit.append(m_n)
			jacobi_sweeps_fit.append(jacobi_run[1])
		else:
			ns_no_fit.append(m_n)
			jacobi_sweeps_no_fit.append(jacobi_run[1])

		sor_sweeps.append(sor_run[1])

		if(m_L <= tipping_point1-L_step0):
			m_L += L_step0
		else:
			m_L += L_step1

	# create the ndarrays to plot
	n_fit_ndarray = np.array(ns_fit)
	jacobi_sweeps_fit_ndarray = np.array(jacobi_sweeps_fit)
	if(jacobi_n_fit_cut != -99):
		n_no_fit_ndarray = np.array(ns_no_fit)
		jacobi_sweeps_no_fit_ndarray = np.array(jacobi_sweeps_no_fit)

	n_ndarray = np.array(ns_fit+ns_no_fit)
	sor_sweeps_ndarray = np.array(sor_sweeps)

	########################################################
	# Plotting

	# Set our colors
	jacobi_color = 'blue'
	sor_color = 'green'

	# Set up the figures and axes
	fig = plt.figure('fig')
	jacobi_ax = fig.add_subplot(111)

	jacobi_ax.set_title('$N_{\mathrm{iter}}(n)$')
	jacobi_ax.set_xlabel('$n = L^{2}$')
	jacobi_ax.set_ylabel('$N_{\mathrm{iter}}$ Jacobi', color=jacobi_color)

	sor_ax = jacobi_ax.twinx() # Get a new axes
	sor_ax.set_ylabel('$N_{\mathrm{iter}}$ SOR', color=sor_color, rotation=-90, labelpad=20)

	# Save handles for legend
	Nn_legend_handles = []

	# Create the plots
	jacobi_data = jacobi_ax.scatter(n_fit_ndarray, jacobi_sweeps_fit_ndarray, marker='o', label='Jacobi $N_{\mathrm{iter}}$', c=jacobi_color)
	Nn_legend_handles.append(jacobi_data)
	if(jacobi_n_fit_cut != -99):
		jacobi_ax.scatter(n_no_fit_ndarray, jacobi_sweeps_no_fit_ndarray, marker='o', c=jacobi_color)
		
	sor_data = sor_ax.scatter(n_ndarray, sor_sweeps_ndarray, marker='d', label='SOR $N_{\mathrm{iter}}$', c=sor_color)
	Nn_legend_handles.append(sor_data)

	# Make the plot log log
	jacobi_ax.set_xscale('log')
	jacobi_ax.set_yscale('log')
	sor_ax.set_xscale('log')
	sor_ax.set_yscale('log')

	# Adjust the axis

	# x1_jacobi_auto,x2_jacobi_auto,y1_jacobi_auto,y2_jacobi_auto = jacobi_ax.axis()
	# jacobi_ax.set_xlim(x_min_scale*x1_jacobi_auto, x_max_scale*x2_jacobi_auto)
	# jacobi_ax.set_ylim(0, y_max_scale*max(jacobi_sweeps_fit+jacobi_sweeps_no_fit))

	x1_sor_auto,x2_sor_auto,y1_sor_auto,y2_sor_auto = sor_ax.axis()
	# sor_ax.set_xlim(x_min_scale*x1_jacobi_auto, x_max_scale*x2_jacobi_auto)
	sor_ax.set_ylim(y1_sor_auto, (10**1.015)*y2_sor_auto)

	# Color the y axis
        for tl in jacobi_ax.get_yticklabels():
                tl.set_color(jacobi_color)
        for tl in sor_ax.get_yticklabels():
                tl.set_color(sor_color)
        for tl in jacobi_ax.get_yminorticklabels():
                tl.set_color(jacobi_color)
        for tl in sor_ax.get_yminorticklabels():
                tl.set_color(sor_color)

	# Fitting 
	########################################################

	########################################################
	# Define the fit function
	def fit_function(n_data, pow_fit, slope_fit, offset_fit):
	        return offset_fit + slope_fit*pow(n_data, pow_fit)
	# end def fit_function

	# Define the fit function2, no offset
	def fit_function2(n_data, pow_fit, slope_fit):
	        return slope_fit*pow(n_data, pow_fit)
	# end def fit_function2

	# actually perform the fits
	# op_par = optimal parameters, covar_matrix has covariance but no errors on plot so it's incorrect...

	# jacobi_p0 = [1.0, 0.0301, 141.45]
	# sor_p0 = [1.0, 1.0, 0.0]

	jacobi_p02 = [1.0, 0.0301]
	sor_p02 = [1.0, 1.0]

	jacobi_fit_status = True
	sor_fit_status = True
	maxfev=m_maxfev = 2000

	try:
		jacobi_op_par, jacobi_covar_matrix = curve_fit(fit_function2, n_fit_ndarray, jacobi_sweeps_fit_ndarray, p0=jacobi_p02, maxfev=m_maxfev)
	except RuntimeError:
		print sys.exc_info()[1]
		print 'jacobi curve_fit failed, continuing...'
		jacobi_fit_status = False 

	try:
		sor_op_par, sor_covar_matrix = curve_fit(fit_function2, n_ndarray, sor_sweeps_ndarray, p0=sor_p02, maxfev=m_maxfev)
	except RuntimeError:
		print sys.exc_info()[1]
		print 'sor curve_fit failed, continuing...'
		sor_fit_status = False 

	# Create fine grain fit x axis
	x1 = min(ns_fit+ns_no_fit)
	x2 = max(ns_fit+ns_no_fit)
	fit_x_ndarray = np.linspace(x1,x2,150)

	# plot the fits
	if(jacobi_fit_status):
		jacobi_fit_line, = jacobi_ax.plot(fit_x_ndarray, fit_function2(fit_x_ndarray, *jacobi_op_par), ls='solid', label='Jacobi Fit', c=jacobi_color)
		Nn_legend_handles.append(jacobi_fit_line)

	if(jacobi_n_fit_cut != -99):
		# Draw vertical line where fit ends
		jacobi_fit_boundary_line = jacobi_ax.axvline(x=jacobi_n_fit_cut, ls = 'dashed', label='Fit Boundary', c='grey')
		Nn_legend_handles.append(jacobi_fit_boundary_line)

	if(sor_fit_status):
		sor_fit_line, = sor_ax.plot(fit_x_ndarray, fit_function2(fit_x_ndarray, *sor_op_par), ls='solid', label='SOR Fit', c=sor_color)
		Nn_legend_handles.append(sor_fit_line)

	# Write out the fit parameters
	# jacobi_fit_text = 'Jacobi Fit Function: $N_{\mathrm{iter}} = c + b n^{a}$' 
	jacobi_fit_text = 'Jacobi Fit Function: $N_{\mathrm{iter}} = b n^{a}$' 
	if(jacobi_fit_status):
		jacobi_fit_text += '\n$a_{\mathrm{Expected}} = 2$\n$a_{\mathrm{Fit}} =$ %.5f' % (jacobi_op_par[0])
		jacobi_fit_text += '\n$b_{\mathrm{Fit}} =$ %.5f' % (jacobi_op_par[1])
		# jacobi_fit_text += '\n$c_{\mathrm{Fit}} =$ %.5f' % (jacobi_op_par[2])
	else:
		jacobi_fit_text += '\nJacobi Fit Failed'

	# sor_fit_text = '\nSOR Fit Function: $N_{\mathrm{iter}} = c + b n^{a}$' 
	sor_fit_text = '\nSOR Fit Function: $N_{\mathrm{iter}} = b n^{a}$' 
	if(sor_fit_status):
		sor_fit_text += '\n$a_{\mathrm{Expected}} = 1$\n$a_{\mathrm{Fit}} =$ %.5f' % (sor_op_par[0])
		sor_fit_text += '\n$b_{\mathrm{Fit}} =$ %.5f' % (sor_op_par[1])
		# sor_fit_text += '\n$c_{\mathrm{Fit}} =$ %.5f' % (sor_op_par[2])
	else:
		sor_fit_text += '\nSOR Fit Failed'

	jacobi_ax.text(0.025, 0.7, jacobi_fit_text+'\n'+sor_fit_text, bbox=dict(edgecolor='black', facecolor='white', fill=True), size='x-small', transform=jacobi_ax.transAxes)

	# Draw the legend
	jacobi_ax.legend(handles=Nn_legend_handles, loc='upper right', bbox_to_anchor=(0.98, 0.98), borderaxespad=0, fontsize='x-small')

	# Annotate
	ann_text = '$R =$ %2.1f [m], $a =$ %1.1f [m]' % (R, a)
	ann_text += '\n$Q/\epsilon_{0} =$ %1.1f [Vm]' % (Q_over_epsilon0)
	ann_text += '\nConv. Criteria: $A <$ %.1E [V]' % (fixed_accuracy)

	jacobi_ax.text(0.679, 0.035, ann_text, bbox=dict(edgecolor='black', facecolor='white', fill=False), size='x-small', transform=jacobi_ax.transAxes)

	# Adjust the minor ticks
	minorFormatter = LogFormatterExponent(labelOnlyBase=False)
	# jacobi_ax.xaxis.set_minor_formatter( minorFormatter )
	jacobi_ax.yaxis.set_minor_formatter( minorFormatter )
	# sor_ax.yaxis.set_minor_formatter( minorFormatter )

	# Print it out
	make_path(m_path)
	fig.savefig(m_path+'/Niter_vs_n.pdf')	

	fig.clf() # Clear fig for reuse

	print 'plot_N_vs_n completed!!!'
# end def for plot_N_vs_n



########################################################
########################################################
########################################################
# Finally, actually run things!

########################################################
########################################################
# Development Runs 

if(False):
	output_path = './output/development/poisson_dipole'
	debugging = False
	
	# target_Dx = 0.1
	# L = int(round(2*R/target_Dx))
	# if(L % 2 == 0): L += 1
	# Dx = float(2*R/L)
	
	
########################################################
########################################################
# Production Runs for paper 

if(True):
	top_output_path = './output/plots_for_paper_final/poisson_dipole'
	debugging = False

        # Part a
        ########################################################
        print '\nPart a:'
        output_path = top_output_path+'/part_a'

	# Set the parameters
	Dx = 0.05 # Set the spacing of the grid points in [m]
	L = int(round(2*R/Dx)) # find the closest L that will work
	if(L % 2 == 0): L += 1 # ensure L is odd

	# Run and print the sim
	m_run = run_sim(L, Dx, 'jacobi', 'epsilon', 0.00001, -99.0, True, '_best_jacobi')
	plot_V('', 'best_jacobi', 60, output_path, m_run)
	plot_Vr_dipole_axis('', 'best_jacobi', output_path, m_run)
        

	# Part b
        ########################################################
        print '\nPart b:'
        output_path = top_output_path+'/part_b'

	plot_N_vs_epsilon(L, Dx, 1.0e-5, 0.5e-5, 15.0e-5, 15.0e-5, 0.001, 50.0e-5, 0.008, output_path)


        ########################################################
        print '\nPart c:'
        output_path = top_output_path+'/part_c'

	plot_N_vs_n(0.001, 35, 1, 50, 20, 238, output_path, 2125)


	# Extra Material
        ########################################################
        print '\nExtra Material:'
        output_path = top_output_path+'/extra_material'

	# do a run that shows the prop boundary
	Dx = 0.1
	L = int(round(2*R/Dx)) # find the closest L that will work
	if(L % 2 == 0): L += 1 # ensure L is odd
	m_run = run_sim(L, Dx, 'jacobi', 'epsilon', 0.02, -99.0, False, '')
	plot_V('', 'prop_illustration', 60, output_path, m_run)


########################################################
print '\n\nDone!\n'


