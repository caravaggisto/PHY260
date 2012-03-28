'''
Created on Feb 7, 2012

@author: res

Assignment: Oscillatory Motion and Chaos

'''

import numpy as np

class ODE_Solver:
    """
    Superclass for numerical methods solving scalar and vector ODEs
    
      du/dt = f(u, t)
    
    Attributes:
    t: array of time values
    u: array of solution values (at time points t)
    k: step number of the most recently computed solution
    f: callable object implementing f(u, t)
    """
    def __init__( self, f ):
        if not callable( f ):
            raise TypeError( 'f is %s, not a function' % type( f ) )
        # For ODE systems, f will often return a list, but
        # arithmetic operations with f in numerical methods
        # require that f is an array. Let self.f be a function
        # that first calls f(u,t) and then ensures that the
        # result is an array of floats.
        self.f = lambda u, t: np.asarray( f( u, t ), float )

    def advance( self ):
        '''
        Advance solution one time step.
        Note: this method is only implemented in subclasses
        for different numerical techniques (e.g. ForwardEuler, 
        EulerCromer, RungeKutta2 4 or 6, etc.
        '''
        raise NotImplementedError

    def set_initial_condition( self, U0 ):
        if isinstance( U0, ( float, int ) ):  # scalar ODE
            self.neq = 1
            U0 = float( U0 )
        else:                            # system of ODEs
            # NOTE: this assumes u is a column vector.
            U0 = np.asarray( U0 )          # (assume U0 is sequence)
            self.neq = U0.size
        self.U0 = U0

    def solve( self, time_points, terminate = None ):
        """
        Compute solution u for t values in the list/array
        time_points, as long as terminate(u,t,step_no) is False. 
        terminate(u,t,step_no) is a user-given function
        returning True or False.
        
        returns tuple (u,t) where
            u = np.array, n x neq 
            t = np.array, n x 1 
        where
            n = # of time points in the time_points array
            neq = order of ODE
        """
        if terminate is None:
            terminate = lambda u, t, step_no: False
            # if not `terminate` not provided, a function that always returns False is used.
            # -> don't stop simulation prematurely for any reason.
        if isinstance( time_points, ( float, int ) ):
            raise TypeError( 'ODE_Solver.solve: time_points is not a sequence' )
            # make sure the time_points input is an array
        self.t = np.asarray( time_points )
        n = self.t.size
        # n = number of time points (NOT the order of the ODE. That would be self.neq)
        if self.neq == 1:  # scalar ODEs
            self.u = np.zeros( n )
        else:              # systems of ODEs
            self.u = np.zeros( ( n, self.neq ) )
        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0
        # Time loop
        for k in range( n - 1 ):
            self.k = k
            self.u[k + 1] = self.advance()
            if terminate( self.u, self.t, self.k + 1 ):
                break  # terminate loop over k
        return self.u, self.t
    
    def interpolate( self, t ):
        '''
        interpolates u[t] for t not in the model's time array.
        '''
        try:
            n = self.t.size
        except TypeError:
            raise TypeError, "Please run the simulation first"
        else:
            if ( n < 2 ):
                raise ValueError, "Not enough simulation steps"
            tmin = self.t[0]
            tmax = self.t[n - 1]
            if t < tmin or t > tmax:
                raise ValueError, \
                  "Requested time is outside the simulated interval"
            for k in range( 0, n - 1 ):
                if t > self.t[k]:
                    break
                    # We want k such that t[k] is the earliest time 
                    # after t, so that 
                    #     self.t[k-1] =< t < self.t[k]
            dt = self.t[k] - self.t[k - 1]
            delta = ( t - self.t[k - 1] ) / dt 
            #     = fraction of dt spanned by t after self.t[k-1]
            u_interpolated = self.u[k - 1] * ( 1.0 - delta ) + self.u[k] * delta
            # = a normalized, weighted sum (equivalently, an average) of u[k-1] and u[k]
            return u_interpolated
    
    
    ''' functions for models of sinusoidally driven ODEs '''
    def get_steady_state_time( self, period ):
        '''
        find the time at which the solution u[t] approaches a steady state.
        How to find? Find t such that:
            - the amplitude stabilizes (has a fractional change < eta, where eta << 1)
                * but, how to determine amplitude? time points aren't necessarily be at peaks.
            - the frequency stabilizes
            - the frequency approaches the driving frequency (this case would be subsumed in the previous condition)
        '''
        #TODO: implement, if necessary :)
        raise NotImplementedError
    
    def sinusoidal( self ):
        #TODO: implement, if necessary. 
        # The idea of this method is just to return True if the data 
        # shows a sinusoidal steady state solution. How to quantify
        # this, not sure. maybe by measuring local extrema in data, 
        # seeing if their magnitude and spacing converge as the time
        # series continues.
        return True
    
    def get_amplitude( self ):
        # return the last local maximum in u(t), assuming sinusoidal data.
        #TODO: find function, maybe in numpy, for local maxima.
        # then, iterate backwards through u(t), starting at end
        if self.sinusoidal():
            k = 1
            while self.u[-k - 1] > self.u[-k]:
                # while downward sloping, working from end of data...
                # increase k.
                k += 1
        return self.u[-k]
    
    def get_phase_shift( self ):
        # External data needed:
        #    driving frequency
        #    phase of driving freq
        # Internal data needed:
        #    if data is sinusoidal,
        #        times at which peaks happen
        #        -> phase of data
        # return (phase of data) - (phase of driving)
        raise NotImplementedError
    
    def get_FWHM( self ):
        # return time gap between half-maxima of steady state part.
        if self.sinusoidal():
            # if data is sinusoidal,
            amplitude = self.get_amplitude()
            half_maximum = amplitude * 0.5
            # extrapolate time points at which this value occurs in the data
        raise NotImplementedError
    
    

class ForwardEuler( ODE_Solver ):
    def advance( self ):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        unew = u[k] + dt * f( u[k], t[k] )
        return unew

class EulerCromer( ODE_Solver ):
    def advance( self ):
        u, f, k, t = self.u, self.f, self.k, self.t
        unew = u[k].copy()
        dt = t[k + 1] - t[k]
        unew[1] = u[k][1] + dt * ( f( u[k], t[k] ) )[1]
        #TODO: possible error. Does this work? Are the indices on these copies working okay?
        ucromer = unew.copy()
        unew[0] = u[k][0] + dt * f( ucromer, t[k] )[0]
        return unew

class RungeKutta2( ODE_Solver ):
    def advance( self ):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        dt2 = dt / 2.0
        K1 = dt * f( u[k], t[k] )
        K2 = dt * f( u[k] + 0.5 * K1, t[k] + dt2 )
        unew = u[k] + K2
        return unew

class RungeKutta4( ODE_Solver ):
    def advance( self ):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        dt2 = dt / 2.0
        K1 = dt * f( u[k], t[k] )
        K2 = dt * f( u[k] + 0.5 * K1, t[k] + dt2 )
        K3 = dt * f( u[k] + 0.5 * K2, t[k] + dt2 )
        K4 = dt * f( u[k] + K3, t[k] + dt )
        unew = u[k] + ( 1 / 6.0 ) * ( K1 + 2 * K2 + 2 * K3 + K4 )
        return unew

import sys, os
# BackwardEuler needs to import function Newton from Newton.py:
try:
    from Newton import Newton
except ImportError:
    pass

class BackwardEuler( ODE_Solver ):
    """Backward Euler solver for scalar ODEs."""
    def __init__( self, f ):
        ODE_Solver.__init__( self, f )
        # Make a sample call to check that f is a scalar function:
        try:
            u = np.array( [1] ); t = 1
            value = f( u, t )
        except IndexError:  # index out of bounds for u
            raise ValueError( 'f(u,t) must return float/int' )

    def advance( self ):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k + 1] - t[k]
        
        def F( w ):
            return w - dt * f( w, t[k + 1] ) - u[k]

        dFdw = Derivative( F )
        w_start = u[k] + dt * f( u[k], t[k] )  # Forward Euler step
        unew, n, F_value = Newton( F, w_start, dFdw, N = 30 )
        if k == 0:
            self.Newton_iter = []
        self.Newton_iter.append( n )
        if n >= 30:
            print "Newton's failed to converge at t=%g "\
                  "(%d iterations)" % ( t[k + 1], n )
        return unew


class Derivative:
    def __init__( self, f, h = 1E-9 ):
        self.f = f
        self.h = float( h )

    def __call__( self, x ):
        f, h = self.f, self.h      # make short forms
        return ( f( x + h ) - f( x - h ) ) / ( 2 * h )


# Tests and verifications

def _f1( u, t ):
    return 0.2 + ( u - _u_solution_f1( t ) ) ** 5

def _u_solution_f1( t ):
    """Exact u(t) corresponding to _f1 above."""
    return 0.2 * t + 3

def _verify( f, exact, U0, T, n ):
    t_points = np.linspace( 0, T, n )
    for method_class in ForwardEuler, RungeKutta4, BackwardEuler:
        try:
            method = method_class( f )
        except:
            continue
        method.set_initial_condition( U0 )
        u, t = method.solve( t_points )
        print method_class.__name__,
        u_exact = np.asarray( exact( t ) ).transpose()
        print 'max error:', ( u_exact - u ).max()
        if method_class is BackwardEuler:
            print 'Backward Euler iterations:', method.Newton_iter

if __name__ == '__main__':
    print 'Exact numerical solution:'
    _verify( _f1, _u_solution_f1, U0 = 3, T = 8, n = 4 )
    print '\nOscillating system:'
    _verify( f = lambda u, t: [u[1], -u[0]],
            exact = lambda t: [np.sin( t ), np.cos( t )],
            U0 = [0, 1], T = 4 * np.pi, n = 20 * 4 )
