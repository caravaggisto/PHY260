'''
Created on Feb 7, 2012

@author: res
'''

from numpy import pi, sin, cos
import numpy as np
#from Vector_3D import V3
#from ODE_Solver import EulerSolver, AboveGround
#from matplotlib.pyplot import plot, text, title, xlabel, ylabel, show, figure, \
#	xlim, ylim, figtext
import math
from ODE_Solver_v3 import EulerCromer, RungeKutta2
from matplotlib.pyplot import *

class Pendulum( object ):
	''' Implementation of different models of pendulum trajectory '''
	def __init__( self, U0, omega_D, alpha_D = 0.2, l = 9.8, gamma = 0.25, linear = True ):
		''' Needed:
			theta_0 (rad)
		Optional (= default):
			omega_0 = 0 (rad s^-1)
			l = 9.8 (m)
			gamma = 0.25 (s^-1)
			alpha_D = 0.2 (rad s^-2) '''
		# Initial conditions
		self.U0 = np.asarray( U0 )
		self.theta_0, self.omega_0 = self.U0[0], self.U0[1]
		## a few constants
		# gravity / restoring force parameters
		self.g, self.l, self.linear = 9.8, l, linear
		# damping force parameters
		self.gamma = gamma
		# driving force parameters
		self.alpha_D, self.omega_D = alpha_D, omega_D
		# an empty array to store model results
		self.model_objects, self.results, self.plot_objects = [], [], []
		
	def get_resonance_frequency( self ):
		''' returns the resonance frequency given an osc driving force '''
		g, l, gamma = self.g, self.l, self.gamma
		return math.sqrt( ( g / l ) - 2 * gamma ** 2 )
	
	def terminate( self, u, t, step_no ):
		return False if u[step_no, 2] >= 0 else True
	
	''' Define ODE'''
	def __call__( self, u, t ):
		'''
		du/dt = f(u,t), where
		u = [ u0 ] = [ theta ] = [ position ]
			[ u1 ]   [ omega ]   [ velocity ]
		f = [ du0/dt   ] = [ u1													]
			[ du1 / dt ]   [ accel_driving(t) + accel_gravity(u0) + accel_damp(u1) ] '''
		theta, omega = u[0], u[1]
		f_out = np.zeros_like( u )
		f_out[0] = omega
		f_out[1] = self.accel_driving( t ) + self.accel_gravity( theta ) + self.accel_damp( omega )
		return f_out
	
	''' Define accelerations acting on the pendulum: '''
	def accel_gravity( self, theta ):
		''' acceleration due to gravity for a pendulum.
		TODO: (optimize) what's the best way to determine 
			the method of a class based on a boolean switch 
			like this? Embedding an if-statement seems 
			computationally inefficient to me. '''
		g, l = self.g, self.l
		if self.linear:
			return -( g / l ) * theta
		else:
			return -( g / l ) * sin( theta )
	
	def accel_damp( self, omega ):
		gamma = self.gamma
		return -2 * gamma * omega
	
	def accel_driving( self, t ):
		alpha_D, omega_D = self.alpha_D, self.omega_D
		return alpha_D * sin( omega_D * t )	
	
	def generate_model( self, method, t_step, t_min, t_max ):
		''' possible methods:
		"EulerCromer", "RungeKutta2", ...
		
		TODO: generate time points over which to model
		TODO: make ODE objects for
			Euler-Cromer
			RK2
		TODO: run the simulation
		TODO: graph it
		TODO: save it to file'''
#		valid_models = { 'EulerCromer'  : EulerCromer( self ), \
#						'RungeKutta2' : RungeKutta2( self ) }
		valid_models = { 'EulerCromer'  : EulerCromer, \
						'RungeKutta2' : RungeKutta2 }

		try:
			method_class = valid_models[method]
			model = method_class()
			model.set_initial_condition( self.U0 )
		except:
			raise ValueError, "method must be one of %s" % valid_models.keys()
		# generate time array
		time_points = np.arange( t_min , t_max + t_step, t_step )
#		solve the ODE system with
		model.solve( time_points )
		return model
	# note: results are accessible via model.u, model.t
	
	def integrate( self, dt, method ):
		'''
			dt = time width ( s )
			method = integration method ( Euler, Runge - Kutta 2nd order, etc )
		TODO: completely gut this crap. This belongs in the ODESolver class.
		'''
		# Run the simulation
		for angle in self.initial_angles_in_radians:
			golf_model = EulerSolver( Fg + Far + Fm, m )
			self.golf_model_objects.append( golf_model )
			self.v0 = V3( self.initial_speed * cos( angle ), \
						  self.initial_speed * sin( angle ), \
						  0.0 )
			golf_model.run( self.x0, self.v0, dt, AboveGround() )
			# Extract the trajectory in the XY plane
			xcoords = [vector.x for vector in golf_model.x]
			ycoords = [vector.y for vector in golf_model.x]
			# Save the data for this model in the local results array
			self.results.append( [ zip( xcoords, ycoords ), dt, angle ] )
			print '%6.2f' % xcoords[-1]
	
	def plot_and_save_models( self, \
							  linewidth = 1.0, \
							  colors = ['r', 'y', 'g', 'b'], \
							  text_size = 10, \
							  text_position = 0.43, \
							  text_position_adjustment = 0.09, \
							  file_directory = "/tmp/", \
							  file_name_base = "PHY260_test", \
							  file_types = ['png'], \
							  show_plot = False, \
							  info_text = '' ):	
		'''
		TODO: MOOOOVE this crap to ODESolver !!! 
		TODO: make this into just a "create figures" method
		TODO: then, make separate methods for 
			( 1 ) viewing and 
			( 2 ) saving to file the figures.
		'''
		# Turn on the grid on the plot
		fig = figure( figsize = ( 7.5, 7.5 ) )
		ax = fig.add_subplot( 111 )
		ax.grid( True )
		# Get data for each model and plot it.
		for golf_ball_model, color, angle in zip( self.golf_model_objects, colors, self.initial_angles_degrees ):
			# gather data from the model object
			xcoords = [vector.x for vector in golf_ball_model.x]
			ycoords = [vector.y for vector in golf_ball_model.x]
			dt = golf_ball_model.dt
			# Draw the trajectory
#			color = colors[ k % len( colors ) ]
#			'''NB: we use a modulus operation in the previous
#			line in case we need to re-use colors. Obviously, this 
#			ideally would be avoided. But, we don't want a silly 
#			error like that stopping the entire process.'''
			plot( xcoords, ycoords, linewidth = linewidth, color = color )
			# Put a label on the plot which shows dt with the same color
			# as the plot line
			#TODO: change the order that these labels appear.
			figtext( .6, text_position + 0.05, "dt = %s, angle = %s" % ( dt, angle ), fontsize = 10, ha = 'left', color = color )
#			text( text_size, text_position, "dt = %s, angle = %s" % ( dt, angle ), fontsize = 10, ha = 'left', color = color )
			text_position += text_position_adjustment
			method = golf_ball_model.name()
		# Put a few finishing touches on the plot
		title( 'Simulated Trajectories of a Golf Ball (%s Method)' % method )
		figtext( .5, .85, info_text, fontsize = 10, ha = 'center' )
		xlabel( 'x (m)' )
		ylabel( 'y (m)' )
		xlim( 0, 500 )
		ylim( 0, 150 )
		fig.patch.set_facecolor( 'white' )
		# Display the results
		if show_plot:
			fig.show()
		for file_type in file_types:
			fig.savefig( file_directory + file_name_base + "." + file_type, facecolor = fig.get_facecolor(), edgecolor = 'none' )

	
