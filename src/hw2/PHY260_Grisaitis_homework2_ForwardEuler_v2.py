'''
Created on Jan 27, 2012

@author: res
'''
import numpy as np

class ForwardEuler_v2( object ):
    """
    Class for solving a scalar of vector ODE,
    
        du/dt = f(u, t)
    
    by the ForwardEuler method.
    
    Class attributes:
    t: array of time values
    u: array of solution values (at time points t)
    k: step number of the most recently computed solution
    f: callable object implementing f(u, t)
    
    Note: This version, v2, of ForwardEuler generalizes v1 to:
        - solve, by Forward Euler, systems of ODEs (not just scalar ODEs)
        - have a user-defined time array (gone with n, T, and dt = T / (n-1) )
        
    """
    def __init__( self, f ):
        if not callable( f ):
            raise TypeError( 'f is % s, not a function' % type( f ) )
        self.f = lambda u, t: np.asarray( f( u, t ) )
    
    def set_initial_condition( self, U0 ):
        if isinstance( U0, ( float, int ) ):    # scalar ODE
            self.neq = 1
        else:                                   # system of ODEs
            # converts the input array to numpy array, if not one already:
            U0 = np.asarray( U0 )
            self.neq = U0.size
        self.U0 = U0
    
    def solve( self, time_points ):
        '''
        Compute u for time points given in time_points.
        
        Note: in this Forward_Euler class, time_points is not given
        at the time of class call. This is in contrast to 
        Forward_Euler_v1, where n and T are given, whence dt (spacing) 
        and the time array t are calculated. 
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
        return self.u, self.t
    
    def advance( self ):
        '''Advance the solution one time step.'''
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        unew = u[k] + dt * f( u[k], t[k] )
        return unew



