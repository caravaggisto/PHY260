'''
Created on Jan 28, 2012

@author: res
'''

import numpy as np

class ODESolver( object ):
    '''
    A superclass to implement multiple numerical techniques for 
    solving ODEs. Instead of having separate classes for each technique, 
    we combine them all here to share the same __init__() 
    (constructor), set_initial_condition() and solve() methods. 
    
    This class combines the functionality that is common to all 
    numerical methods in approximating solutions to ODEs (and systems
    thereof), such as:
        * u = numpy array holding the solutions at discrete time points
        * t = numpy array of time points
        * f = a callable function object that holds information about the ODE
        * k = index number corresponding to the current time step
        * U0 = initial condition for u
        * a method to implement the time-step approximations for 
        each numerical technique.
            - this is done with 2 methods, actually:
                1. solve() for performing the time loop
                2. advance() for advancing the time steps
            - the advance() method is not defined in the superclass.
            - advance() is the bit that varies by technique.
        * 
    '''
    
    def __init__( self, f ):
        # f is the function such that dN/dt = f(N, t)
        if not callable( f ):
            raise TypeError( 'f is % s, not a function' % type( f ) )
        self.f = lambda u, t: np.asarray( f( u, t ), float )
    
    def advance( self ):
        '''
        NOTE: this only does the Forward Euler approximation 
        technique. I'm saving implementation of other numerical
        techniques for another day, if it's necessary.
        '''
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        unew = u[k] + dt * f( u[k], t[k] )
        return unew
    
    '''
#    This code is only to be used when we want to employ different
#    numerical techniques. For homework #2, we're only using the 
#    Forward Euler algorithm.
    def advance( self ):
        # This method is only implemented in subclasses corresponding to
        # different numerical solving methods.
        raise NotImplementedError
    '''
    
    def set_initial_condition( self, U0 ):
        if isinstance( U0, ( float, int ) ):    # scalar ODE
            self.neq = 1
            U0 = float( U0 )
        else:                                   # system of ODEs
            # converts the input array to numpy array, if 
            # it's not one already
            U0 = np.asarray( U0 )
            self.neq = U0.size
        self.U0 = U0
    
    def solve( self, time_points, terminate = None ):
        '''
        Compute u for time points given in time_points.
        
        Note #1: in this ODESolver class, time_points is not given
        at the time of class call. This is in contrast to 
        Forward_Euler_v1, where n and T are given, whence dt ( spacing ) 
        and the time array t are calculated. 
        
        Note #2: the terminate parameter is a function that returns a
        boolean, telling the solver whether or not to continue 
        approximating. For example, if we're simulating a ball 
        falling, we don't care about predictions when it's 
        "underground". terminate() can prevent such wasteful 
        computations from occurring. 
        '''
        # First, store the time points into the self object. This way, 
        # time steps will be visible in the advance() method (see below).
        self.t = np.asarray( time_points )
        # n, however, is not needed anywhere else, so we make it local:
        n = self.t.size
        if self.neq == 1:                       # scalar ODE
            # u will store approximations for u(t)
            self.u = np.zeros( n )
        else:                                   # system of ODEs
            # u will store approx's for u(t), u'(t), u''(t), ...
            self.u = np.zeros( ( n, self.neq ) )
        # Assume self.t[0] corresponds to self.U0
        self.u[0] = self.U0
        # Time loop
        for k in range( n - 1 ):
            self.k = k
            self.u[k + 1] = self.advance()
            if terminate is not None:
                if terminate( self.u, self.t, self.k + 1 ):
                    break   # terminate loop after k.
        return self.t, self.u
    
    def plot_results( self ):
#        TODO: implement! Look in firefox bookmarks for Python - ODE Computation - Zombie Apocalypse example
        pass




#'''
#sub-classes for implementing different numerical techniques,
#other than Forward Euler...
#'''
#
#class Forward_Euler( ODESolver ):
#    def advance( self ):
#        u, f, k, t = self.u, self.f, self.k, self.t
#        dt = t[k + 1] - t[k]
#        unew = u[k] + dt * f( u[k], t[k] )
#        return unew
