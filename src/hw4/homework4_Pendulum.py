'''
Created on Feb 7, 2012

@author: res
'''

from numpy import sin, sqrt
import numpy as np
import math

class Pendulum( object ):
	''' Implementation of the ODE describing a pendulum: '''
#	def __init__( self, omega_D, alpha_D = 0.2, l = 9.8, gamma = 0.25, linear = True ):
	def __init__( self, omega_D, alpha_D, l, gamma, linear ):
		''' Parameters:
			omega_D = driving frquency		(rad s-1)
			l 		= length of pendulum	(m)
			gamma 	= damping coefficient	(s-1)
			alpha_D = driving amplitude		(rad s-2) '''
		# gravity / restoring force parameters
		self.g, self.l, self.linear = 9.8, l, linear
		# damping force parameters
		self.gamma = gamma
		# driving force parameters
		self.alpha_D, self.omega_D = alpha_D, omega_D
		
	def get_resonance_frequency( self ):
		''' returns the resonance frequency given an osc driving force '''
		g, l, gamma = self.g, self.l, self.gamma
		return sqrt( ( g / l ) - 2 * gamma ** 2 )
	
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
