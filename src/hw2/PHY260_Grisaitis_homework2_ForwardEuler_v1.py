'''
Created on Jan 27, 2012

@author: res
'''
import numpy as np

class ForwardEuler_v1( object ):
    def __init__( self, f, U0, T, n ):
        '''
        Constructor
            f = a callable function; the first derivative of U
            U0 = initial condition for U
            T = time cap
            n = number of steps, after initial value.
            dt = time gap per step
            u = array of U(t)'s
            t = array of time points.
        '''
        self.f, self.U0, self.T, self.n = f, U0, T, n
        self.dt = T / float( n )
        self.u = np.zeros( n + 1 )
        self.t = np.zeros( n + 1 )
    
    def solve( self ):
        """Compute solution for 0 <= t <= T."""
        self.u[0] = float( self.U0 )
        self.t[0] = float( 0 )
        for k in range( self.n ):
            self.k = k
            self.t[k + 1] = self.t[k] + self.dt
            self.u[k + 1] = self.advance()
        return self.u, self.t
    
    def advance( self ):
        """Advance the solution one time step.
        Note how this isolates the numerical process. That is, 
        this is where the Euler-ness is happening.
        """
        u, dt, f, k, t = self.u, self.dt, self.f, self.k, self.t
        unew = u[k] + dt * f( u[k], t[k] )
        return unew


