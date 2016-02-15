#!/usr/bin/env /home/mbe9/Documents/Spring_2016/PHYS_566_Computational_HW/anaconda/bin/python
# #! to the correct python install, not portable!!!

# Run by calling
# ./chaotic_pendulum.py 2>&1 | tee output.log
# to save output to output.log

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

########################################################
# Set fixed parameters [SI units]
g = 9.80665 # g, gravity [m/s^2]
l = 9.8 # l, length [m]
gamma = 2.05 # gamma, dampening constant [s^-1]
possible_alphaD = [0.2, 0.5, 1.2] # driving force amplitude [rad/s^2]

# Euler-Cromer method parameters
Dt = 0.001 # Time step of TODO simulation
max_Dtheta = np.pi/100.0 # TODO Maximum radial distance allowed between time steps

# TODO list of 10 good driving frequencies for res plot

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
    # There are no real private variables in python so I can't enforce this...
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
	run[0].tmax = periodsD*(2.0*np.pi/omegaD)
	run[0].sim_method = sim_method
	run[0].lin_case = lin_case
	
	# Set the initial values
	run[0].set_theta(theta0)
	run[0].omega = 0.0

	# Loop for periodsD driving force periods
	n = 0 # time step index


	while( n*Dt <= run[0].tmax ):
		# append a new time_step object for n+1
		run.append(time_step((n+1)*Dt)) 

		if(sim_method == 'euler_cromer'): # Euler-Cromer method

			# Compute the new velocity, see writeup for details...
			run[n+1].omega = run[n].omega + ( -(g/l)*linearity(run[n].theta, lin_case) -2*gamma*run[n].omega + alphaD*np.sin(omegaD*run[n].t) )*Dt
	
			# Compute the new position
			run[n+1].set_theta( run[n].theta + run[n+1].omega*Dt )

		#if(sim_method == 'runge_kutta'): # Runge-Kutta method

			# TODO



		# Make sure we didn't move to far, notify user and exit if we did
		Dtheta = abs(run[n+1].rawtheta - run[n].theta) # Use rawtheta to avoid branch cut problems
		if( Dtheta > max_Dtheta):
			print 'On time step %d the pendulum moved to far, Dtheta = %.3f\nProgram Exiting' % (n+1, Dtheta )
			sys.exit()

		n = n+1 # increment time step
	# end while loop

	print 'Run Completed!'

	return run
# end def for run_sim



########################################################
# Define a function to take two runs and compare them
# Runs should have the same driving force and run time
# Should be used to compare Euler-Cromer to Runge-Kutta and linear to nonlinear
def compare_runs(run1, run1_name, run2, run2_name, title, fname):
	if( run1[0].alphaD != run2[0].alphaD and run1[0].omegaD != run2[0].omegaD and run1[0].tmax != run2[0].tmax ):
		print 'SERIOUS ERROR: Comparing incompatible runs!!!!'

	########################################################
	# Create ndarrays that we can plot
	num_timesteps = len(run1)
	
	t_ndarray = np.zeros(num_timesteps)

	driving_theta = np.zeros(num_timesteps)

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

		driving_theta[i] = run1[0].alphaD*np.sin(run1[0].omegaD*run1[i].t) 

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

	plt.plot(t_ndarray, run1_theta, ls='solid', label=run1_name, c='blue')
	plt.plot(t_ndarray, run2_theta, ls='solid', label=run2_name, c='green')
	plt.plot(t_ndarray, driving_theta, ls='dotted', label='$F_{Driving}$', c='black')

	theta_ax.set_title(title)
	theta_ax.set_xlabel('$t$ [s]')
	theta_ax.set_ylabel('$\\theta$ [radians]')
	theta_ax.set_xlim((0.0, run1[0].tmax))
	plt.legend()

	theta_fig.savefig(m_path+'/'+fname+'_theta.pdf')

	# run1 energy
	run1_energy_fig = plt.figure('run1_energy') # get a separate figure
	run1_energy_ax = run1_energy_fig.add_subplot(111) # Get the axes, effectively

	plt.plot(t_ndarray, run1_kenetic, ls='solid', label='$K$', c='blue')
	plt.plot(t_ndarray, run1_potential, ls='solid', label='$U$', c='green')
	plt.plot(t_ndarray, run1_total, ls='solid', label='$E$', c='black')

	run1_energy_ax.set_title('Energy: '+run1_name)
	run1_energy_ax.set_xlabel('$t$ [s]')
	run1_energy_ax.set_ylabel('$E$ [J/mass]')
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
	run2_energy_ax.set_ylabel('$E$ [J/mass]')
	run2_energy_ax.set_xlim((0.0, run1[0].tmax))
	plt.legend()

	run2_energy_fig.savefig(m_path+'/'+fname+'_'+run2_name+'_energy.pdf')




	print 'compare_runs completed'
# end def for compare_runs



########################################################
########################################################
########################################################
# Finally, actually run things!

output_path = './output' # Set output path, hard coded...

# vary theta0
'''
# run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case)

ec_run = run_sim(possible_alphaD[0], 1.0, 0.0, 25, 'euler_cromer','linear')
ec_run2 = run_sim(possible_alphaD[0], 1.0, 20.0*(np.pi/180.0), 25, 'euler_cromer','linear')

# compare_runs(run1, run1_name, run2, run2_name, title, fname) 
compare_runs(ec_run, '$\\theta_{0} = 0.0^{\circ}$', ec_run2, '$\\theta_{0} = 10.0^{\circ}$', 'Vary $\\theta_{0}$', 'vary_theta0')
'''

# vary linearity 
# run_sim(alphaD, omegaD, theta0, periodsD, sim_method, lin_case)

ec_run = run_sim(possible_alphaD[0], 1.0, np.pi, 25, 'euler_cromer','linear')
ec_run2 = run_sim(possible_alphaD[0], 1.0, np.pi, 25, 'euler_cromer','nonlinear')

# compare_runs(run1, run1_name, run2, run2_name, title, fname) 
compare_runs(ec_run, 'Linear', ec_run2, 'Nonlinear', 'Vary Linearity', 'vary_linearity')




########################################################
print '\n\nDone!\n'


