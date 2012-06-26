'''
Created on Jan 27, 2012

@author: res
'''
import numpy as np
from numpy import log as ln
from numpy import exp
from PHY260_Grisaitis_homework2_ODESolver import ODESolver
#from matplotlib.pyplot import plot
#import matplotlib.pyplot as plt


class Radioactive_Decay( object ):
    ''' Mathematical model for the radioactive decay ODE:
        dN/dt = - (   1   /    tau    ) * N
              = - ( ln(2) / half_life ) * N
    Attributes to initialize:
        initial_mass    in kg
        molar_mass      in g per mol
        half_life       in years
    Attributes made later:
        N0              initial moles, derived from initial_mass and molar_mass
        N               array of approximations for N(t)
        time_points     array of time points used in last approximation 
    Methods:
        __call__(N,t)   the callable function f(N,t) where dN/dt = f(N,t)
        __str__()       print f(N,t) in readable text
        approximation_Forward_Euler( time_array )
                        returns ( N(t), t ) numerical approximation given a 
                        time_array of points at which to evaluate.
        reset_initial_conditions( initial_mass )
                        ... self-explanatory
    
    Notes:
        Making this function as a class, and not a stand-alone function,
        lets us do operations on it that only pertain to certain para-
        meters (e.g. differentiate or integrating w.r.t. time, change 
        half-life, etc.)
    '''
    
    def initial_moles( self ):
        # Number of moles: (mass in kg) * ( 1000g/1kg ) / (g / mol)
        return float( self.initial_mass * 1000 / self.molar_mass )
    
    def __init__( self, initial_mass, molar_mass, half_life ):
        ''' initialize the object:
        initial_mass    in kg
        molar_mass      in g per mol
        half_life       in years
        '''
        self.initial_mass = initial_mass
        self.molar_mass = molar_mass
        self.half_life = float( half_life )
        self.N0 = self.initial_moles()
        # make some empty arrays for results data.
#        self.N = np.array( () )
        self.solutions = []
#        self.time_points = np.array( () )
        self.time_arrays = []
#        figure_number keeps track of when we create new figures for graphing.
#        self.figures = []
#        self.figure_number = 0
    
    def __call__( self, N, t ):
        ''' return dN/dt for the radioactive decay ODE'''
        half_life = self.half_life
        # Returns a scalar, because this is a system of one ODE.
        return -N * ln( 2 ) / half_life
    
    def reset_initial_conditions( self, initial_mass ):
        self.initial_mass = initial_mass
        self.set_initial_moles()
    
    def approximation_Forward_Euler( self, time_points ):
        
        '''
        approximate N(T) via the Forward Euler technique
        time_points = array of time points at which to approximate.
        returns N, t
        '''
        self.time_arrays.append( time_points )
        # initialize an ODE solver to approximate our system.
        forward_euler_approximation = ODESolver( self.__call__ )
        forward_euler_approximation.set_initial_condition( self.N0 )
        t, N = forward_euler_approximation.solve( time_points )
        # store the result in the object
        self.solutions.append( N )
        return t, N
    
    def analytical_solution( self, t ):
        ''' Returns the exact solution N(t), given
        N0              initial condition, scalar
        t               point in time, scalar 
        '''
        return self.N0 * exp( -ln( 2 ) * t / self.half_life )
    
    def analytical_solution_array( self, time_points ):
        '''
        returns (t, N(t)), an array of *exact* solution points given
        time_points, a numpy array of time points at which to evaluate 
        N(t).
        '''
        num_data_points = self.time_points.size
        exact_results = np.zeros( num_data_points )
        for k in range( num_data_points ):
            exact_results[k] = self.analytical_solution( time_points[k] )
        return time_points, exact_results
    
    def deviation_from_exact( self, after_time_T, which_solution ):
#        TODO: fill in!!!
        pass
    
