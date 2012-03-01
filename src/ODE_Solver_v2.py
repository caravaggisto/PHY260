'''
Author: William C Grisaitis
'''

"""
This module implements several algorithms for solving particle motion ODEs
"""

import math
import numpy as np

class OdeSolver_v2( object ):
    """
    Base class for particle motion ODE solving schemes. The subclasses
    must implement the following function:

    _step(self, dt, t, x, v) should return a tuple dx, dv which are
    the coordinate and velocity increments for evolving from the current
    step to the next.

    Also, the subclasses must define the class variable _name which
    should be set to a descriptive string naming the method.

    The subclasses may also provide the function _setup(self) which can
    be used for additional initialization.
    """
    def __init__( self, f ):
        '''
        args:
            f(u,t)
                u is a np array 
                t is a scalar
                returns np array
        '''
        self.f = f              # f needs to be a callable function like f(u,t) that returns an ARRAY
        self.u_initial = None   # u needs to be an array, where du/dt = f(u,t)
        self.u = None
        self.t = None
        self._setup()
    
    def _setup( self ):
        pass
    
    def name( cls ):
        return cls._name
    name = classmethod( name )

    # Pythonic convention for implementing "pure virtual"
    # methods is to raise the "NotImplementedError" exception
    def _step( self, dt, t, x, v ):
        raise NotImplementedError, "Function not implemented"

    def run( self, u_initial, dt, runningCondition ):
        """
        Run the simulation and accumulate the history.
        Arguments are as follows:

        xInitial         -- Initial coordinates
        vInitial         -- Initial velocities
        dt               -- Simulation time step
        runningCondition -- Callable which should return "true" at each
                            step if we are to continue the simulation.
                            It is called as runningCondition(t, x, v)
        """
        # Make sure dt is float. Then t will be float as well.
        dt = dt * 1.0
        self.dt = dt
        # Running the simulation backwards requires special care
        # which is beyond our purposes here
        assert dt > 0.0, "Can not run the simulation backward in time"
        # Remember the initial conditions and initialize the history
#        self.x0 = xInitial * 1.0
#        self.v0 = vInitial * 1.0
        self.u_initial = u_initial
        #TODO: make these numpy arrays. Then, change whatever methods need to be changed when they're used in `run` below
        self.t = np.asarray( [0.0] )
#        self.t = [0.0, ]
        self.u = np.asarray( [u_initial] )
#        self.u = [u_initial, ]
        # Initialize the running variables
        nsteps = 0
        t = 0.0
        u = self.u_initial
        while ( runningCondition( t, u ) ):
            #TODO: make running_condition take argument `u` (np array), NOT x, v
#            dx, dv = self._step( dt, t, x, v )
            du = self._step( dt, t, u )
#            x = x + dx
#            v = v + dv
            u = u + du
            # If we were not accumulating the history,
            # the use of += operator would be more appropriate,
            # for example x += dx. However, the += operator
            # does not necessarily create a new object, and we
            # do not want to fill the history with references
            # to the same object! This illustrates, again,
            # the distinction between Python variables and
            # variables in other programming languages.
            #
            # Also, we do not use t = t + dt because this leads
            # to accumulation of round-off errors. For x and v
            # we have no choice but for t we can avoid it.
            nsteps += 1
            t = dt * nsteps
            self.t = np.append( self.t, [t] )
#            self.t.append( t )
            #TODO: when self.u is made into a np array, change these `appends`
            self.u.append( u )
#            self.x.append( x )
#            self.v.append( v )   
         
    
    def evolve( self, u_initial, dt, runningCondition ):
        """
        Run the simulation without accumulating the history.
        Arguments are as follows:
        
        u_initial        -- (np array) Initial coordinates, velocities, etc
        dt               -- (float) Simulation time step
        runningCondition -- (Callable) which should return "true" at each
                            step if we are to continue the simulation.
                            It is called as runningCondition(t, u)
        """
        # Make sure dt is float. Then t will be float as well.
        dt = dt * 1.0
        # Running the simulation backwards requires special care
        # which is beyond our purposes here
        assert dt > 0.0, "Can not run the simulation backward in time"
        # Remember the initial conditions
        self.u_initial = u_initial * 1.0
        # Initialize the running variables
        nsteps = 0
        t = 0.0
        u = self.u_initial
        while ( runningCondition( t, u ) ):
            du = self._step( dt, t, u )
#            dx, dv = self._step( dt, t, x, v )
            u += du
#            x += dx
#            v += dv
            nsteps += 1
            t = dt * nsteps
        self.t = t
        self.u = u    
    
    def interpolate( self, t ):
        """
        This function can be invoked after calling "run" in order to
        interpolate system coordinates and velocity to an arbitrary time
        moment covered by the simulation. The interpolation is linear
        between the two history points closest to the given time. Note that
        linear interpolation can significantly degrade the resolution of
        a high-order method such as 4th order Runge-Kutta. You should
        not use this function to sample the simulation history at time
        intervals shorter than the simulation time step -- instead, just
        rerun the simulation using a smaller step.
        """
        try:
            n = len( self.t )
        except TypeError:
            # self.t is not a sequence. Re-raise the exception
            # with an appropriate error message.
            raise TypeError, "Please run the simulation first"
        else:
            if ( n < 2 ):
                raise ValueError, "Not enough simulation steps"
            tmin = self.t[0]
            tmax = self.t[n - 1]
            if t < tmin or t > tmax:
                raise ValueError, \
                      "Requested time is outside the simulated interval"
            dt = ( tmax - tmin ) * 1.0 / ( n - 1 )
            nbelow = int( math.floor( ( t - tmin ) / dt ) )
            nabove = nbelow + 1
            if nabove >= n:
                nabove = n - 1
                nbelow = nabove - 1
            delta = ( t - tmin ) / dt - nbelow
            u = self.u[nbelow]
            #TODO: change to use arrays
            x = self.x[nbelow] * ( 1.0 - delta ) + self.x[nabove] * delta
            v = self.v[nbelow] * ( 1.0 - delta ) + self.v[nabove] * delta
            return x, v

###########################################################################
#
# Some concrete implementations of the base class follow
#
###########################################################################

class EulerSolver( OdeSolver_v2 ):
    """
    Implements the Euler ODE solving scheme. Construct the solver with
    EulerSolver(F, m) where F is the force model and m is the particle mass.
    The force model is a callable which will be called by the solver like
    this: F(t, x, v).

    It is also possible to call the solver as EulerSolver(A) where A is
    the acceleration as a function of t, x, and v (unit mass is assumed).
    """
    _name = "Euler"
#    def _step( self, dt, t, x, v ):
#        # Euler scheme: x and v are taken at the current t
#        dv = dt / self.mass * self.Force( t, x, v )
#        dx = v * dt
#        return dx, dv
    def _step( self, dt, t, u ):
        # Euler scheme: x and v are taken at the current t
        du = dt * self.f( u, t )
        return du

class EulerCromer( OdeSolver_v2 ):
    """
    Implements the Euler-Cromer ODE solving scheme. Construct the solver with
    EulerCromer(F, m) where F is the force model and m is the particle mass.
    The force model is a callable which will be called by the solver like
    this: F(t, x, v).

    It is also possible to call the solver as EulerCromer(A) where A is
    the acceleration as a function of t, x, and v (unit mass is assumed).
    """
    _name = "Euler-Cromer"
#    def _step( self, dt, t, x, v ):
#        # Note the difference with the Euler method: the value of dx
#        # is calculated using the new value of v, not the old one
#        self.dt = dt
#        dv = dt / self.mass * self.Force( t, x, v )
#        dx = ( v + dv ) * dt
#        return dx, dv
    def _step( self, dt, t, u ):
        # Note the difference with the Euler method: the value of dx
        # is calculated using the new value of v, not the old one
        self.dt = dt
        #TODO: beware of array([1,2]) or array([[1,2]])
        du = np.array( [0, 0] )
        du[1] = dt * self.f( u, t )[1]
        du[0] = dt * ( u[1] + du[1] )
        return du

class RK2( OdeSolver_v2 ):
    '''
    TODO: implement!!!
    The 2nd-order Runge-Kutta method:
    u[k+1] = u[k] + K2
    K1 = dt * f(u[k], t[k])
    K2 = dt * f(u[k] + 0.5*K1, t[k] + 0.5*dt) 
    '''
    #TODO: update to work with arrays, u, self.f(u,t)
    def _setup( self ):
        ''' note how this self._f is basically making my self.f '''
        # For simplicity, define a function which combines
        # x and v variables into one tuple
#        self._f = lambda t, u: ( v[], self.f( t, u ) )
    
    def _step( self, dt, t, u ):
        f = self.f
#        halfstep = dt / 2.0
#        k1x, k1v = self._f( t, x, v )
#        k2x, k2v = self._f( t + halfstep, x + halfstep * k1x, v + halfstep * k1v )
#        k3x, k3v = self._f( t + halfstep, x + halfstep * k2x, v + halfstep * k2v )
#        k4x, k4v = self._f( t + dt, x + dt * k3x, v + dt * k3v )
#        dx = dt / 6.0 * ( k1x + 2 * k2x + 2 * k3x + k4x )
#        dv = dt / 6.0 * ( k1v + 2 * k2v + 2 * k3v + k4v )
        return du


class RK4( OdeSolver_v2 ):
    """
    Implements the 4th order Runge-Kutta ODE solving scheme. Construct the
    solver with RK4(F, m) where F is the force model and m is the particle
    mass. The force model is a callable which will be called by the solver
    like this: F(t, x, v).

    It is also possible to call the solver as RK4(A) where A is the
    acceleration as a function of t, x, and v (unit mass is assumed).
    """
    _name = "4th order Runge-Kutta"
    #TODO: update to work with arrays, u, self.f(u,t)
    def _setup( self ):
        # For simplicity, define a function which combines
        # x and v variables into one tuple
        self._f = lambda t, x, v: ( v, self.Force( t, x, v ) / self.mass )
    
    def _step( self, dt, t, x, v ):
        halfstep = dt / 2.0
        k1x, k1v = self._f( t, x, v )
        k2x, k2v = self._f( t + halfstep, x + halfstep * k1x, v + halfstep * k1v )
        k3x, k3v = self._f( t + halfstep, x + halfstep * k2x, v + halfstep * k2v )
        k4x, k4v = self._f( t + dt, x + dt * k3x, v + dt * k3v )
        dx = dt / 6.0 * ( k1x + 2 * k2x + 2 * k3x + k4x )
        dv = dt / 6.0 * ( k1v + 2 * k2v + 2 * k3v + k4v )
        return dx, dv


class RK6( OdeSolver_v2 ):
    """
    Implements the 6th order Runge-Kutta ODE solving scheme. Construct the
    solver with RK6(F, m) where F is the force model and m is the particle
    mass. The force model is a callable which will be called by the solver
    like this: F(t, x, v).

    It is also possible to call the solver as RK6(A) where A is the
    acceleration as a function of t, x, and v (unit mass is assumed).
    """
    _name = "6th order Runge-Kutta"
    #TODO: update to work with arrays, u, self.f(u,t)

    def _setup( self ):
        # For simplicity, define a function which combines
        # x and v variables into one tuple
        self._f = lambda h, t, x, v: ( h * v, h / self.mass * self.Force( t, x, v ) )
    
    def _step( self, dt, t, x, v ):
        h = dt * 1.0
        kx1, kv1 = self._f( h, t, x, v )
        kx2, kv2 = self._f( h, t + h / 3, x + kx1 / 3, v + kv1 / 3 )
        kx3, kv3 = self._f( h, t + 2 * h / 3, x + 2 / 3.0 * kx2, v + 2 / 3.0 * kv2 )
        kx4, kv4 = self._f( h, t + h / 3, x + ( kx1 + 4 * kx2 - kx3 ) / 12.0,
                                     v + ( kv1 + 4 * kv2 - kv3 ) / 12.0 )
        kx5, kv5 = self._f( h, t + h / 2, x + ( 18 * kx2 - 3 * kx3 - 6 * kx4 - kx1 ) / 16.0,
                                     v + ( 18 * kv2 - 3 * kv3 - 6 * kv4 - kv1 ) / 16.0 )
        kx6, kv6 = self._f( h, t + h / 2, x + ( 9 * kx2 - 3 * kx3 - 6 * kx4 + 4 * kx5 ) / 8.0,
                                     v + ( 9 * kv2 - 3 * kv3 - 6 * kv4 + 4 * kv5 ) / 8.0 )
        kx7, kv7 = self._f( h, t + h, x + ( 9 * kx1 - 36 * kx2 + 63 * kx3 + 72 * kx4 - 64 * kx5 ) / 44.0,
                                   v + ( 9 * kv1 - 36 * kv2 + 63 * kv3 + 72 * kv4 - 64 * kv5 ) / 44.0 )
        dx = ( 11 * ( kx1 + kx7 ) + 81 * ( kx3 + kx4 ) - 32 * ( kx5 + kx6 ) ) / 120.0
        dv = ( 11 * ( kv1 + kv7 ) + 81 * ( kv3 + kv4 ) - 32 * ( kv5 + kv6 ) ) / 120.0
        return dx, dv


###########################################################################
#
# Some common running/stopping conditions for the ODE solvers
#
###########################################################################

class TimeLimit:
    def __init__( self, tmax ):
        self.tmax = tmax
    def __call__( self, t, u ):
        #TODO: change the orders of these arguments
        return t < self.tmax

class AboveGround:
    def __call__( self, t, u ):
        #TODO: change the orders of these arguments
        return u[0] >= 0.0
