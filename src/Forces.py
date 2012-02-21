'''
This code mostly comes from http://highenergy.phys.ttu.edu/, where 
Igor Volobouev provided these classes to handle physical forces for 
use in a Computational Physics course at Texas Tech.
'''

from Vector_3D import V3
import math
from math import exp

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
        sum1 = None
        for k, val in self._basisCallables.iteritems():
            if sum1 is None:
                sum1 = k( t, x, v ) * val
            else:
                sum1 += k( t, x, v ) * val
        return sum1

########################################################################
#
# Concrete implementations of various forces follow
#
########################################################################

class Force_Gravity( BasicForce ):
    def __init__( self, g, m ):
        self.force_vector = V3( 0.0, -g * m, 0.0 )
    def __call__( self, t, x, v ):
        return self.force_vector

class Magnus_Effect( BasicForce ):
    def __init__( self, S0, omega_vector ):
        self.S0_omega = S0 * omega_vector
    def __call__( self, t, x, v ):
        return self.S0_omega.cross( v )

class Ball_Drag( BasicForce ):
    '''
    For the drag, assume a general form of
        Fdrag = -C rho A v**2
    '''
    def __init__( self, A, rho, C, dimpled ):
        self.A = float( A )
        self.rho = float( rho )
        self.C = float( C )
        self.dimpled = dimpled
    
    def __call__( self, t, x, v ):
        speed = v.length()
        if speed > 14.0 and self.dimpled:
            C = 7.0 / speed
        else:
            C = self.C
        return -C * self.A * self.rho * speed * v
