"""
This module implements Euler's algorithm for solving particle motion ODEs.
"""

__author__="Igor Volobouev (i.volobouev@ttu.edu)"
__version__="0.1"
__date__ ="Jan 22 2008"

class EulerSolver:
    """
    Implements the Euler ODE solving scheme. Construct the solver with
    EulerSolver(F, m) where F is the force model and m is the particle mass.
    The force model is a callable which will be called by the solver like
    this: F(t, x, v).

    It is also possible to call the solver as EulerSolver(A) where A is
    the acceleration as a function of t, x, and v (unit mass is assumed).
    """
    def __init__(self, F, m=1.0):
        self.Force = F
        self.mass = m*1.0
        if (self.mass <= 0.0):
            raise NotImplementedError, \
                  "Photons and exotic matter are not supported"
        self.x0 = None
        self.v0 = None
        self.t = None
        self.x = None
        self.v = None

    def run(self, xInitial, vInitial, dt, runningCondition):
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
        self.x0 = xInitial*1.0
        self.v0 = vInitial*1.0
        t = 0.0
        x = self.x0
        v = self.v0
        self.t = [t,]
        self.x = [x,]
        self.v = [v,]
        while (runningCondition(t, x, v)):
            # Euler scheme: x and v are taken at the current t
            dv = self.Force(t, x, v)/self.mass*dt
            x = x + v*dt
            v = v + dv
            t = t + dt
            # If we were not accumulating the history,
            # the use of += operator would be more appropriate,
            # for example x += v*dt. However, the += operator
            # does not necessarily create a new object, and we
            # do not want to fill the history with references
            # to the same object! This illustrates, again,
            # the distinction between Python variables and
            # variables in other programming languages.
            self.t.append(t)
            self.x.append(x)
            self.v.append(v)

    def evolve(self, xInitial, vInitial, dt, runningCondition):
        """
        Run the simulation without accumulating the history.
        Arguments are as follows:

        xInitial         -- Initial coordinates
        vInitial         -- Initial velocities
        dt               -- Simulation time step
        runningCondition -- Callable which should return "true" at each
                            step if we are to continue the simulation.
                            It is called as runningCondition(t, x, v)
        """
        self.x0 = xInitial*1.0
        self.v0 = vInitial*1.0
        t = 0.0
        x = self.x0
        v = self.v0
        while (runningCondition(t, x, v)):
            dv = self.Force(t, x, v)/self.mass*dt
            x += v*dt
            v += dv
            t += dt
        self.t = t
        self.x = x
        self.v = v
