#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

# Run by calling
# ./chaotic_pendulum.py 2>&1 | tee output.log
# to save output to output.log

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

########################################################
# Set fixed parameters [SI units]
g = 9.8 # g, gravity [m/s^2]
l = 9.8 # l, length [m]
gamma = 0.25 # gamma, dampening constant [s^-1]
possible_alphaD = [0.2, 0.5, 1.2] # driving force amplitude [rad/s^2]

# Euler-Cromer method parameters
Dt = 0.001 # Time step of TODO simulation
max_Dtheta = np.pi/100.0 # TODO Maximum radial distance allowed between time steps

# compute omega0
omega0 = g/l

# compute resonant omega from theory 
omega_res_theory = np.sqrt( omega0 - 2*np.square(gamma) )


########################################################
# Print out starting values
print '\nBeginning chaotic_pendulum.py simulation'
print '\nFixed Parameters are:'
print '---------------------------------------------'
print 'g = %.5f [m/s^2]' % g
print 'l = %.5f [m]' % l
print 'Dampening Constant gamma = %.5f [s^-1]' % gamma

print '\nPossible Driving Amplitudes alphaD [radians/s^2]'
print possible_alphaD

print '\nomega0 = %.2f [radians/s]' % omega0
print 'omega res theory = %.5f [radians/s]' % omega_res_theory

print '\nSimulation Time Step = %.5f [s]' % Dt
print 'Simulation Max Step Size = %.3f [radians]' % max_Dtheta
print '---------------------------------------------'
print '---------------------------------------------\n'


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
# Define time_step class to hold all the relevant parameters at a time step
class time_step:

    def __init__(self, time):
	self.t = time
        self.theta = -99.0
	self.omega = -99.0
	self.forceD = -99.0

	# save the raw theta coordinate, no bounds
	self.rawtheta = -99.0

	# save the sim parameters, only fill at t = 0
	self.alphaD = -99.0
	self.omegaD = -99.0
	self.theta0 = -99.0
	self.tmax = -99.0
	self.sim_method = ''
	self.lin_case = ''

    # Define a function to set theta while forcing -pi < theta <= pi 
    # WARNING do not use object.theta = X to set theta!!!
    # I may clean this up using python properties, at the very end... TODO 
    def set_theta(self, new_theta):
	self.rawtheta = new_theta

	while ( new_theta <= -np.pi) : new_theta = new_theta + 2.0*np.pi
	while ( new_theta > np.pi) : new_theta = new_theta - 2.0*np.pi

	self.theta = new_theta
	return None # explicitly return none, python would implicitly otherwise

# end class for time_step

########################################################
# Define a function to return theta or sin(theta) depending on the lin_case 
def linearity(theta, lin_case):
	if(lin_case == 'linear'):
		return theta
	elif(lin_case == 'nonlinear'):
		return np.sin(theta)
# end def for linearity

########################################################
# Define a function to run the simulation
def run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case):
	
	# Start the list of time_steps
	run = [time_step(0.0)]

	# save the sim parameters
	run[0].alphaD = alphaD
	run[0].omegaD = omegaD
	run[0].theta0 = theta0
	run[0].tmax = periodsD*((2.0*np.pi)/omegaD)
	run[0].sim_method = sim_method
	run[0].lin_case = lin_case
	
	# Set the initial values
	run[0].set_theta(theta0)
	run[0].omega = 0.0
	run[0].forceD = alphaD*np.sin(omegaD*run[0].t)

	# Loop for periodsD driving force periods
	n = 0 # time step index


	while( n*Dt <= run[0].tmax ):
		# append a new time_step object for n+1
		run.append(time_step((n+1)*Dt)) 

		run[n+1].forceD = alphaD*np.sin(omegaD*run[n+1].t)

		if(sim_method == 'euler_cromer'): # Euler-Cromer method

			# Compute the new velocity, see writeup for details...
			run[n+1].omega = run[n].omega + ( -omega0*linearity(run[n].theta, lin_case) -2*gamma*run[n].omega + alphaD*np.sin(omegaD*run[n].t) )*Dt
	
			# Compute the new position
			run[n+1].set_theta( run[n].theta + run[n+1].omega*Dt )

		#if(sim_method == 'runge_kutta'): # Runge-Kutta method

			# TODO



		# Make sure we didn't move to far, notify user and exit if we did
		Dtheta = abs(run[n+1].rawtheta - run[n].theta) # Use rawtheta to avoid branch cut problems
		if( Dtheta > max_Dtheta):
			print 'On time step %d the pendulum moved to far, Dtheta = %.3f\nProgram Exiting' % (n+1, Dtheta )
			sys.exit()

		n += 1 # increment time step
	# end while loop

	print 'Run Completed!'

	return run
# end def for run_sim



########################################################
# Define a function to take two runs and compare them
# Runs should have the same driving force and run time
# Should be used to compare Euler-Cromer to Runge-Kutta and linear to nonlinear, or different starting angles
def compare_runs(run1, run1_name, run2, run2_name, title, fname):
	print 'Beginning Compare Runs:'
	if( run1[0].alphaD != run2[0].alphaD and run1[0].omegaD != run2[0].omegaD and run1[0].tmax != run2[0].tmax ):
		print 'SERIOUS ERROR: Comparing incompatible runs!!!!'

	########################################################
	# Create ndarrays that we can plot
	num_timesteps = len(run1)
	
	t_ndarray = np.zeros(num_timesteps)

	driving_force = np.zeros(num_timesteps)

	run1_theta = np.zeros(num_timesteps)
	run1_omega = np.zeros(num_timesteps)
	run1_kenetic = np.zeros(num_timesteps)
	run1_potential = np.zeros(num_timesteps)
	run1_total = np.zeros(num_timesteps)

	run2_theta = np.zeros(num_timesteps)
	run2_omega = np.zeros(num_timesteps)
	run2_kenetic = np.zeros(num_timesteps)
	run2_potential = np.zeros(num_timesteps)
	run2_total = np.zeros(num_timesteps)

	# Fill the ndarrays
	mass = 1.0 # put mass in formulas just out of habit
	# energy units will be energy/mass
	for i in range(num_timesteps):
		t_ndarray[i] = run1[i].t

		driving_force[i] = run1[i].forceD

		run1_theta[i] = run1[i].theta
		run1_omega[i] = run1[i].omega
		run1_kenetic[i] = 0.5*mass*l*l*np.square(run1[i].omega)
		run1_potential[i] = mass*g*l*(1.0-np.cos(run1[i].theta))
		run1_total[i] = run1_kenetic[i] + run1_potential[i]

		run2_theta[i] = run2[i].theta
		run2_omega[i] = run2[i].omega
		run2_kenetic[i] = 0.5*mass*l*l*np.square(run2[i].omega)
		run2_potential[i] = mass*g*l*(1.0-np.cos(run2[i].theta))
		run2_total[i] = run2_kenetic[i] + run2_potential[i]


	########################################################
	# Create the plots
	m_path = output_path+'/'+fname
	make_path(m_path)

	# theta
	theta_fig = plt.figure('theta') # get a separate figure
	theta_ax = theta_fig.add_subplot(111) # Get the axes, effectively

	theta_ax.set_title(title)
	theta_ax.set_xlabel('$t$ [s]')
	theta_ax.set_ylabel('$\\theta$ [radians]')

	theta_ax.set_xlim((0.0, run1[0].tmax))

	force_theta_ax = theta_ax.twinx() # Get a new axes for the force
	force_theta_color = 'darkmagenta'
	force_theta_ax.set_ylabel('$F_{Driving}$', color=force_theta_color, rotation=-90)
	force_theta_ax.set_ylim(-1.5*run1[0].alphaD, 1.5*run1[0].alphaD)

	theta_ax.plot(t_ndarray, run1_theta, ls='solid', label=run1_name, c='blue')
	theta_ax.plot(t_ndarray, run2_theta, ls='solid', label=run2_name, c='green')
	force_theta_ax.plot(t_ndarray, driving_force, ls='dotted', label='$F_{Driving}$', c=force_theta_color)

	x1_theta,x2_theta,y1_theta,y2_theta = theta_ax.axis()
	theta_ax.set_ylim(-1.2*max(y1_theta,y2_theta), 1.2*max(y1_theta,y2_theta))

	theta_ax.legend()

	for tl in force_theta_ax.get_yticklabels():
		tl.set_color(force_theta_color)

	theta_fig.savefig(m_path+'/'+fname+'_theta.pdf')

	# omega
	omega_fig = plt.figure('omega') # get a separate figure
	omega_ax = omega_fig.add_subplot(111) # Get the axes, effectively

	omega_ax.set_title(title)
	omega_ax.set_xlabel('$t$ [s]')
	omega_ax.set_ylabel('$\\omega$ [radians/s]')
	omega_ax.set_xlim((0.0, run1[0].tmax))

	force_omega_ax = omega_ax.twinx() # Get a new axes for the force
	force_omega_color = 'darkmagenta'
	force_omega_ax.set_ylabel('$F_{Driving}$', color=force_omega_color, rotation=-90)
	force_omega_ax.set_ylim(-1.5*run1[0].alphaD, 1.5*run1[0].alphaD)

	omega_ax.plot(t_ndarray, run1_omega, ls='solid', label=run1_name, c='blue')
	omega_ax.plot(t_ndarray, run2_omega, ls='solid', label=run2_name, c='green')
	force_omega_ax.plot(t_ndarray, driving_force, ls='dotted', label='$F_{Driving}$', c=force_omega_color)

	x1_omega,x2_omega,y1_omega,y2_omega = omega_ax.axis()
	omega_ax.set_ylim(-1.2*max(y1_omega,y2_omega), 1.2*max(y1_omega,y2_omega))

	omega_ax.legend()

	for tl in force_omega_ax.get_yticklabels():
		tl.set_color(force_omega_color)

	omega_fig.savefig(m_path+'/'+fname+'_omega.pdf')


	# run1 energy
	run1_energy_fig = plt.figure('run1_energy') # get a separate figure
	run1_energy_ax = run1_energy_fig.add_subplot(111) # Get the axes, effectively

	plt.plot(t_ndarray, run1_kenetic, ls='solid', label='$K$', c='blue')
	plt.plot(t_ndarray, run1_potential, ls='solid', label='$U$', c='green')
	plt.plot(t_ndarray, run1_total, ls='solid', label='$E$', c='black')

	run1_energy_ax.set_title('Energy: '+run1_name)
	run1_energy_ax.set_xlabel('$t$ [s]')
	run1_energy_ax.set_ylabel('$E$/mass [J/mass]')
	run1_energy_ax.set_xlim((0.0, run1[0].tmax))
	plt.legend()

	run1_energy_fig.savefig(m_path+'/'+fname+'_'+run1_name+'_energy.pdf')


	# run2 energy
	run2_energy_fig = plt.figure('run2_energy') # get a separate figure
	run2_energy_ax = run2_energy_fig.add_subplot(111) # Get the axes, effectively

	plt.plot(t_ndarray, run2_kenetic, ls='solid', label='$K$', c='blue')
	plt.plot(t_ndarray, run2_potential, ls='solid', label='$U$', c='green')
	plt.plot(t_ndarray, run2_total, ls='solid', label='$E$', c='black')

	run2_energy_ax.set_title('Energy: '+run2_name)
	run2_energy_ax.set_xlabel('$t$ [s]')
	run2_energy_ax.set_ylabel('$E$/mass [J/mass]')
	run2_energy_ax.set_xlim((0.0, run1[0].tmax))
	plt.legend()

	run2_energy_fig.savefig(m_path+'/'+fname+'_'+run2_name+'_energy.pdf')




	print 'Compare Runs completed'
# end def for compare_runs

########################################################
# Define a function to run and fit one omegaD for resonance sweeping
def res_run(omegaD, alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case, m_path):

	# Perform a run 
	run = run_sim(alphaD, omegaD, theta0, run_end_periodsD, sim_method, lin_case)

	# set up fit ranges
	fit_range_min = fit_begin_periodsD*(2.0*np.pi/omegaD)
	fit_range_max = run[0].tmax
	
	########################################################
	# As we can't specify a fit range when we fit, we have to force it by
	# dividing up our data into chunks by hand. ROOT can do this much better,
	# scipy failed me! First we create python lists as they can be built up from
	# nothing, then we make ndarrays so numpy works as intended. It should be able
	# to plot and fit lists just fine, but it was having weird behavior... 
	# It's quite a hack but if it works it works...

	# Declare the lists
	t_driving_force = []
	driving_force = []
	t_nofit = []
	theta_nofit = []
	t_fit = []
	theta_fit = []

	# Fill the lists
	for i in range(len(run)):
		t_driving_force.append(run[i].t)
		driving_force.append(run[i].forceD)

		if( run[i].t < fit_range_min or run[i].t > fit_range_max):
			t_nofit.append(run[i].t)
			theta_nofit.append(run[i].theta)
		elif( run[i].t >= fit_range_min and run[i].t <= fit_range_max):
			t_fit.append(run[i].t)
			theta_fit.append(run[i].theta)
			# Fill nofit with Nones so it doesn't plot at the y axis min
			# we aren't fitting this one so it's fine
                        t_nofit.append(None)
                        theta_nofit.append(None)
	# Create the ndarrays

	t_driving_force_ndarray = np.array(t_driving_force)
	driving_force_ndarray = np.array(driving_force)
	t_nofit_ndarray = np.array(t_nofit)
	theta_nofit_ndarray = np.array(theta_nofit)
	t_fit_ndarray = np.array(t_fit)
	theta_fit_ndarray = np.array(theta_fit)

	########################################################
	# Create the theta plot
	fig = plt.figure('res_run') # get a separate figure
	theta_ax = fig.add_subplot(111) # Get the axes, effectively

	theta_title = '$\\theta\left(t\\right)$, $\Omega_{D} = $ %.4f [radians/s] = %.2f $\Omega_{Res Theory}$' % (omegaD, omegaD/omega_res_theory) 
	theta_ax.set_title(theta_title)
	theta_ax.set_xlabel('$t$ [s]')
	theta_ax.set_ylabel('$\\theta$ [radians]')

	# Get a new axes for the force so we can plot it with a different y scale
	force_theta_ax = theta_ax.twinx()
 	force_theta_color = 'darkmagenta'
	force_theta_ax.set_ylabel('$F_{Driving}$', color=force_theta_color, rotation=-90)
	# we known alphaD should be the force y max, so set y axis range now
	force_yax_multiplier = 1.6 # control the force y axis range
	force_theta_ax.set_ylim(-force_yax_multiplier*alphaD, force_yax_multiplier*alphaD)
	for tl in force_theta_ax.get_yticklabels():
		tl.set_color(force_theta_color)
	
	# Create the three base plots
	theta_ax.plot(t_nofit_ndarray, theta_nofit_ndarray, ls='dotted', label='Unfitted $\\theta\left(t\\right)$', c='blue')
	theta_ax.plot(t_fit_ndarray, theta_fit_ndarray, ls='solid', label='Fitted $\\theta\left(t\\right)$', c='blue')
	force_theta_ax.plot(t_driving_force_ndarray, driving_force_ndarray, ls='dotted', label='$F_{Driving}$', c=force_theta_color)

	# Now that they have been plotted, clean up the theta y range
	x1_theta,x2_theta,y1_theta,y2_theta = theta_ax.axis()
	theta_yax_multiplier = 1.5 # control the theta y axis range
	theta_ax.set_ylim(-theta_yax_multiplier*max(y1_theta,y2_theta), theta_yax_multiplier*max(y1_theta,y2_theta))

	# Fitting 
	########################################################
	# Define a sinusoidal fit function
	def sine_fit_function(theta, thetaP_fit, omegaD_fit, phiP_fit):
		return thetaP_fit*np.sin(omegaD_fit*theta - phiP_fit)
	# end def sine_fit_function

	# Set fit color
	fit_color = 'green'

	# Add vertical line at where the fit range begins
	theta_ax.axvline(x=fit_range_min, ls = 'solid', label='Fit Begins', c='gray')

	# compute particular solution parameters from theory
	thetaP_theory = alphaD/np.sqrt(np.square(np.square(omega0) - np.square(omegaD)) + 4*np.square(gamma)*np.square(omegaD) ) 
	phiP_theory = np.arctan((2*gamma*omegaD)/(np.square(omega0) - np.square(omegaD)))
	# omegaP theory is just omegaD...

	# set them as the initial/guess fit parameters
	m_p0 = [thetaP_theory, omegaD, phiP_theory]

	# actually perform the fit
	# op_par = optimal parameters, covar_matrix has covariance but no errors on plot so it's incorrect...
	op_par, covar_matrix = curve_fit(sine_fit_function, t_fit, theta_fit, p0=m_p0)

	# get the number of data points, to set the marker spacing by doing integer division // on it
	num_points = len(t_fit)

	# plot the fit
	theta_ax.plot(t_fit_ndarray, sine_fit_function(t_fit_ndarray, *op_par), ls='None', markevery=num_points//70, marker='s', markersize=7, label='Fit', c=fit_color)

	# plot the theory particular/steady state solution
	theta_ax.plot(t_fit_ndarray, sine_fit_function(t_fit_ndarray, *m_p0), ls='None', markevery=num_points//70, marker='d', markersize=6, label='$\\theta\left(t\\right)_{Particular}$', c='maroon')

	# Write out the fit parameters
	fit_text = '$\\theta_{P Theory} =$ %.5f, $\\theta_{P Fit} =$  %.5f' % (m_p0[0], op_par[0])
	fit_text += '\n$\\Omega_{D} =$ %.5f, $\\Omega_{D Fit} =$ %.5f' % (m_p0[1], op_par[1])
	fit_text += '\n$\phi_{Theory} =$ %.5f, $\phi_{Fit} =$ %.5f' % (m_p0[2], op_par[2])
	plt.figtext(0.61, 0.13, fit_text, bbox=dict(edgecolor='black', fill=False), size='x-small' )

	# draw final legend, set the x range, and print it out!
	theta_ax.legend(bbox_to_anchor=(0.025, 0.92, 0.925, 0.10), loc=3, ncol=5, mode="expand", borderaxespad=0.0, fontsize='small')

	# full view
	theta_ax.set_xlim((0.0, run[0].tmax))
	fname = 'res_run_omegaD_over_omega_res_theory_%.2f.pdf' % (omegaD/omega_res_theory)
	fig.savefig(m_path+'/'+fname)

	# fit region only
	theta_ax.set_xlim((0.9*fit_range_min, 1.05*fit_range_max))
	fname = 'res_run_omegaD_over_omega_res_theory_%.2f_fit_region.pdf' % (omegaD/omega_res_theory)
	fig.savefig(m_path+'/'+fname)

	# Clear the figure for the next res_sweep
	fig.clf()

	print 'Resonance Run Completed!'

	# Return the list of fit parameters
	return op_par

# end def for res_run

########################################################
# Define a function to sweep for resonances
def res_sweep(num_runs, omegaD_percent_range, alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case):
	print 'Beginning Resonance Sweep:'

	# make paths
	m_path = output_path+'/res_sweep'
	make_path(m_path)

	m_path2 = m_path+'/res_runs'
	make_path(m_path2)

	# make num runs odd so we always hit omega res theory in the middle
	if(num_runs % 2 == 0): num_runs = num_runs + 1

	# make a ndarray of omegaD's to sweep over, our spectrum
	m_omegaD_min = (1.0-omegaD_percent_range)*omega_res_theory
	m_omegaD_max = (1.0+omegaD_percent_range)*omega_res_theory
	omegaD_spectrum_ndarray = np.linspace(m_omegaD_min, m_omegaD_max, num=num_runs)

	# Create empty list to store lists of fit parameters
	all_op_pars = []

	# Create empty list to store fit amplitudes, op_par[0]
	thetaP_fits = []

	# Create empty list to store fit phases, op_par[2]
	phi_fits = []

	# run res_run for each omegaD
	for i in range(omegaD_spectrum_ndarray.size):
		all_op_pars.append(res_run(omegaD_spectrum_ndarray[i], alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case, m_path2))
		thetaP_fits.append(all_op_pars[i][0])
		phi_fits.append(all_op_pars[i][2])

	# Note you must manually view each res_run output plot to make sure the fit isn't bad
	# TODO could implement auto Rsquared check of some kind... very end

	# Create ndarrays of the amplitudes and phases we can plot nicely, no bugs!
	thetaP_fits_ndarray = np.array(thetaP_fits)
	phi_fits_ndarray = np.array(phi_fits)


	########################################################
	# Define a plotting function for our spectrum sweeps 
	def spectrum_plot_function(y_variable_ndarray, y_variable_name, fname):

		fig = plt.figure('res_sweep_fig') # get a separate figure
		ax = fig.add_subplot(111) # Get the axes, effectively

		ax.set_title(y_variable_name+' vs $\Omega_{D}$')
		ax.set_xlabel('$\Omega_{D}$ [radians/s]')
		ax.set_ylabel(y_variable_name+' [radians]')

		# Create the base plot
		ax.scatter(omegaD_spectrum_ndarray, y_variable_ndarray, marker='o', label=y_variable_name, c='blue')

		# Add a vertical line at omega res theory
		ax.axvline(x=omega_res_theory, ls = 'dashed', label='$\Omega_{Resonance Theory}$', c='gray')


		'''
		# Fitting 
		########################################################
		# Define a TODO fit function
		def sine_fit_function(theta, thetaP_fit, omegaD_fit, phiP_fit):
			return thetaP_fit*np.sin(omegaD_fit*theta - phiP_fit)
		# end def sine_fit_function

		# Set fit color
		fit_color = 'green'
	
		# set them as the initial/guess fit parameters
		m_p0 = [thetaP_theory, omegaD, phiP_theory]

		# actually perform the fit
		# op_par = optimal parameters, covar_matrix has covariance but no errors on plot so it's incorrect...
		op_par, covar_matrix = curve_fit(sine_fit_function, t_fit, theta_fit, p0=m_p0)

		# plot the fit
		ax.plot(t_fit_ndarray, sine_fit_function(t_fit_ndarray, *op_par), ls='None', markevery=num_points//70, marker='s', markersize=7, label='Fit', c=fit_color)

		# Write out the fit parameters
		fit_text = '$\\theta_{P Theory} =$ %.5f, $\\theta_{P Fit} =$  %.5f' % (m_p0[0], op_par[0])
		plt.figtext(0.61, 0.13, fit_text, bbox=dict(edgecolor='black', fill=False), size='x-small' )
		'''
		# Write out the simulation parameters
		sim_text = '$\\alpha_{D} =$ %.1f [radians/$s^2$], $\\theta_{0} =$  %.1f$^{\circ}$' % (alphaD, theta0*(180.0/np.pi))
		if( sim_method == 'euler_cromer'): m_sim_method = 'Euler-Cromer'
		if( sim_method == 'runge_kutta'): m_sim_method = 'Runge-Kutta'
		if( lin_case == 'linear'): m_lin_case = 'Linear'
		if( lin_case == 'nonlinear'): m_lin_case = 'Nonlinear'
		sim_text += '\nSim. Method = %s\nLinearity = %s' % (m_sim_method, m_lin_case)
		plt.figtext(0.15, 0.13, sim_text, bbox=dict(edgecolor='black', fill=False), size='x-small' )

		# draw final legend
		# ax.legend(bbox_to_anchor=(0.025, 0.92, 0.925, 0.10), loc=3, ncol=5, mode="expand", borderaxespad=0.0, fontsize='small')
		ax.legend()

		# set the x range
		#ax.set_xlim((0.0, run[0].tmax))

		# print it out
		fig.savefig(m_path+'/res_sweep_'+fname+'.pdf')

		# Clear the figure for the next spectrum
		fig.clf()

		# end def spectrum_plot_function 

	########################################################

	# Make amplitude vs omegaD plot
	spectrum_plot_function(thetaP_fits_ndarray, '$\\theta_{P Fit}$', 'thetaP')

	# Make phase vs omegaD plot
	spectrum_plot_function(phi_fits_ndarray, '$\phi_{Fit}$', 'phi')

	print 'Resonance Sweep Completed!!!!'

# end def for res_sweep

########################################################
########################################################
########################################################
# Finally, actually run things!

output_path = './output'

# vary theta0
'''
# run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case)

ec_run = run_sim(possible_alphaD[0], 1.0, 0.0, 25, 'euler_cromer','linear')
ec_run2 = run_sim(possible_alphaD[0], 1.0, 20.0*(np.pi/180.0), 25, 'euler_cromer','linear')

# compare_runs(run1, run1_name, run2, run2_name, title, fname) 
compare_runs(ec_run, '$\\theta_{0} = 0.0^{\circ}$', ec_run2, '$\\theta_{0} = 10.0^{\circ}$', 'Vary $\\theta_{0}$', 'vary_theta0')
'''

# vary linearity 
'''
# run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case)

ec_run = run_sim(possible_alphaD[0], 1.0, np.pi, 25, 'euler_cromer','linear')
ec_run2 = run_sim(possible_alphaD[0], 1.0, np.pi, 25, 'euler_cromer','nonlinear')

# compare_runs(run1, run1_name, run2, run2_name, title, fname) 
compare_runs(ec_run, 'Linear', ec_run2, 'Nonlinear', 'Vary Linearity', 'vary_linearity')
'''

# raw res_runs
'''
# All of these need to be given a m_path if run again...
# res_run(omegaD, alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case)
res_run(0.8*np.sqrt( omega0 - 2*np.square(gamma) ), possible_alphaD[0], 0.0, 6, 10, 'euler_cromer', 'linear')
res_run(1.1*np.sqrt( omega0 - 2*np.square(gamma) ), possible_alphaD[0], 0.0, 6, 10, 'euler_cromer', 'linear')
res_run(1.5*np.sqrt( omega0 - 2*np.square(gamma) ), possible_alphaD[0], 0.0, 6, 10, 'euler_cromer', 'linear')
'''

# res_sweep
# res_sweep(num_runs, omegaD_percent_range, alphaD, theta0, fit_begin_periodsD, run_end_periodsD, sim_method, lin_case)

# res_sweep(3, 0.6, possible_alphaD[0], 0.0, 12, 16, 'euler_cromer', 'linear')

res_sweep(33, 0.8, possible_alphaD[0], 0.0, 12, 16, 'euler_cromer', 'linear')








########################################################
print '\n\nDone!\n'


