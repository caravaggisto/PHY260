'''
Created on Apr 15, 2012

@author: willg
'''

import numpy as np
from datetime import datetime
from time import time

def montecarlo_integration_v1( f, a, b, n ):
    ''' monte carlo integration of
        f = callable function
        a = start coordinate (an int)
        b = end coordinate
        n = number of iterations for the approximation (resolution of the MC integral)
    TODO: make this work for vector functions.
    '''
    s = 0 # = total for the integration, to be accumulated through multiple iterations.
    for i in range( n ):
        x = np.random.uniform( a, b ) # random number in [a,b).
        s += f( x )
    I = float( b - a ) * s / n
    return I # = float, approximated integral of f on [a,b).

def montecarlo_integration_v2( f, a, b, n ):
    ''' same as v1, but with vectorized computation (no `for` loop)
    --> ~10x faster
    '''
    x = np.random.uniform( a, b, n )
    s = np.sum( f( x ) )
    I = float( b - a ) * s / n
    return I # = float, approximated integral of f on [a,b).

def montecarlo_integration_v3( f, a, b, n, N = 100 ):
    ''' monte carlo integration of
        f = callable function
        a = start coordinate (an int)
        b = end coordinate
        n = number of iterations for the approximation (resolution of the MC integral)
    TODO: make this work for vector functions.
    '''
    s = 0 # = total for the integration, to be accumulated through multiple iterations.
    I_values = np.zeros( int( n / N ) ) # = array of results of I.
    k_values = np.zeros_like( I_values )
    for k in range( 1, n + 1 ):
        x = np.random.uniform( a, b )
        s += f( x ) # = running total of integrations
        if k % N == 0:
            I_values[( k / N ) - 1] = float( b - a ) * s / k # = average of the first k integrations.
            k_values[( k / N ) - 1] = k
    return k_values, I_values

class test_func( object ):
    ''' a callable math function '''
    def __call__( self, x ):
        return 2 + 3 * x

def func( x ):
    return 2 + 3 * x

if __name__ == '__main__':
    a, b = 1, 2
    n = int( 1e7 )
    for MC_integrator in ( montecarlo_integration_v1, montecarlo_integration_v2 ):
        time1 = time()
        print "%.10f in %.5f seconds" % ( MC_integrator( func, a, b, n ), time() - time1 )
