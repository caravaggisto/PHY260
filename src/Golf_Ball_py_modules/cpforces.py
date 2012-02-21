"""
This module implements physical force models which can be used with ODE
solvers in the Computational Physics course.
"""

__author__ = "Igor Volobouev (i.volobouev@ttu.edu)"
__version__ = "0.2"
__date__ = "Jan 29 2008"

import v3
import math
from math import exp

# Forces should behave like vectors. There is a complication however:
# apriori, we do not know the dimensionality of the velocity space
# and we do not know how many times we need to calculate the forces.
#
# This problem can be solved efficiently with the following design:
# have all force classes inherit from the "BasicForce" class which
# implements the vector behavior. Each "BasicForce" instance works as
# a basis vector in the vector space of forces at all simulation times.
# The dimensionality of this space is D*N where D is the dimensionality of
# the velocity space (usually 3) and N is the number of simulation steps.
#
# The "CompositeForce" class which we will create as a result of
# arithmetic operations with forces will represent a linear superposition
# of basis vectors. In this way we will need to run the calculations
# for each force only once at each simulation step, no matter how many
# times it is used in a linear combination.
#
class BasicForce( object ):
    "Base class for all forces in ODE solvers."
    # Derived classes should override the __call__ function
    def __call__( self, t, x, v ):
        raise NotImplementedError, "Operation not implemented"
    # __hash__ and __cmp__ functions will allow us to use
    # BasicForce instances as dictionary keys
    def __hash__( self ):
        return id( self )
    def __cmp__( self, other ):
        return id( self ) - id( other )
    # Subsequent functions implement basic vector operations
    def __add__( self, other ):
        return CompositeForce( 1.0, self, 1.0, other )
    def __sub__( self, other ):
        return CompositeForce( 1.0, self, -1.0, other )
    def __mul__( self, other ):
        return CompositeForce( other, self, 0.0, self )
    def __rmul__( self, other ):
        return CompositeForce( other, self, 0.0, self )
    def __div__( self, other ):
        return self * ( 1.0 / other )
    def __neg__( self ):
        return CompositeForce( -1.0, self, 0.0, self )
    def __pos__( self ):
        return CompositeForce( 1.0, self, 0.0, self )

class CompositeForce( BasicForce ):
    def __init__( self, c1, f1, c2, f2 ):
        self._basisCallables = dict()
        if ( c1 != 0.0 or ( c1 == 0.0 and c2 == 0.0 ) ):
            self._include( c1, f1 )
        if ( c2 != 0.0 ):
            self._include( c2, f2 )
    def _include( self, c, f ):
        # Is f itself a composite force?
        if ( hasattr( f, "_basisCallables" ) ):
            for k, v in f._basisCallables.iteritems():
                self._basisCallables[k] = self._basisCallables.get( k, 0.0 ) + c * v
        else:
            self._basisCallables[f] = self._basisCallables.get( f, 0.0 ) + c
    def __call__( self, t, x, v ):
        sum = None
        for k, val in self._basisCallables.iteritems():
            if sum is None:
                sum = k( t, x, v ) * val
            else:
                sum += k( t, x, v ) * val
        return sum

########################################################################
#
# Concrete implementations of various forces follow
#
########################################################################

class ForceOfGravity( BasicForce ):
    "The force of gravity near the Earth surface"
    def __init__( self, m ):
        self.f = v3.V3( 0.0, -9.8 * m, 0.0 )
    def __call__( self, t, x, v ):
        return self.f

class MagnusForce( BasicForce ):
    "Magnus force assuming rotation with constant agular velocity"
    def __init__( self, S0, omega ):
        self.somega = S0 * omega
    def __call__( self, t, x, v ):
        return self.somega.cross( v )

class QuadraticDrag( BasicForce ):
    "Quadratic drag with constant drag coefficient"
    def __init__( self, dragCoefficient, frontalArea, airDensity ):
        self.C = dragCoefficient
        self.A = frontalArea
        self.rho = airDensity
    def __call__( self, t, x, v ):
        # Return it in a form which works both for scalar
        # and vector velocity. Also note that v is the last
        # factor. This means the formula multiplies scalars
        # first and, when v is a vector, it multiplies a scalar
        # by a vector only at the very end which saves
        # a little bit of time.
        return -0.5 * self.C * self.A * self.rho * abs( v ) * v

class GolfBallDrag( BasicForce ):
    """
    Air drag acting on a golf ball using the model of Giordano
    and Nakanishi (Section 2.5)
    """
    def __init__( self ):
        # Golf ball has a diameter of 42.67 mm
        self.A = math.pi * 42.67e-3 ** 2 / 4.0
        self.rho = 1.2
    def __call__( self, t, x, v ):
        s = abs( v )
        if ( s < 14.0 ):
            C = 0.5
        else:
            C = 7.0 / s
        return -C * self.A * self.rho * s * v

class AltitudeDrag( BasicForce ):
    """
    Quadratic drag with air density dependence on altitude using
    isothermal atmosphere
    """
    def __init__( self, dragCoefficient, frontalArea, seaLevelAirDensity ):
        self.C = dragCoefficient
        self.A = frontalArea
        self.rho0 = seaLevelAirDensity
        self.y0 = 1.0e4
    def __call__( self, t, x, v ):
        rho = self.rho0 * exp( -x.y / self.y0 )
        return -0.5 * self.C * self.A * rho * abs( v ) * v

class RestoringForce( BasicForce ):
    "Restoring force proportional to linear displacement"
    def __init__( self, springConstant, x0 ):
        self.springConstant = springConstant
        self.x0 = x0
    def __call__( self, t, x, v ):
        return -self.springConstant * ( x - self.x0 )
