'''
Created on Apr 22, 2012

@author: willg
'''
import numpy as np

# For debugging, to use same random set.
#np.random.seed( 1 )

class Diffusion( object ):
    ''' Wrapper for the 1D diffusion equation, which we want to...
        - solve via difference methods
    Diffusion is given by:
        du/dt = D*d2u/dx2
            = (D/dx^2)*(u[i+1,n] -2u[i,n] +u[i-1,n)
    or, equivalently,
        u[i,n+1] = u[i,n] + (D*dt/dx^2)*(u[i+1,n] -2u[i,n] +u[i-1,n)
    where
        x = i * dx
        t = n * dt
    '''
    def __init__( self, X, T, dx, dt, D ):
        self.X, self.T, self.dx, self.dt, self.D = X, T, dx, dt, D
        # x's range from 0 to X-1. Same for t's
        self.u = np.zeros( ( X, T ) )
        # and set the initial condition:
        midpoint = X / 2
        for x in [midpoint + offset for offset in [-2, -1, 0, 1]]:
            self.u[x, 0] = 0.25
    
    def advance( self, i, n ):
        ''' returns u[i,n+1] where x = i*dx, t = n*dt, for 
        u[i,n+1] = u[i,n] + (D*dt/dx^2)*(u[i+1,n] -2u[i,n] +u[i-1,n)
        '''
        u, D, dx, dt = self.u, self.D, self.dx, self.dt
        # u[i, n + 1] = ... 
        return u[i, n] + ( D * dt / dx ** 2 ) * ( u[i + 1, n] - 2 * u[i, n] + u[i - 1, n ] )
    
    def solve( self ):
        self.i, self.n = 0, 0 # initial x and t values
        X, T = self.X, self.T
        #TODO: to optimize this, vectorize it. For X, T < 100, 100, it probably doesn't matter.
        for n in range( 0, T - 1 ):
            for i in range( 1, X - 1 ):
                self.u[i, n + 1] = self.advance( i, n )
    '''
    TODO: make methods that compute metadata, stored locally
    '''
    def get_x_mean( self ):
        ''' computes and returns <x> '''
        # note, expect dimensions to be (T,1)
        self.x_means = self.u.mean( axis = 0 )
    
    def get_x2_mean( self ):
        ''' computes and returns <x^2> '''
        self.x2_means = ( self.u ** 2 ).mean( axis = 0 )

if __name__ == '__main__':
    '''this code is for debugging the class, making sure 
    the methods are all working the way they're supposed to, etc. '''
    np.random.seed( 1 ) # use the same set of random numbers... forever!
    #===========================================================================
    X = 100 # number of 1D sites
    T = 1000 # number of time steps
    dx = 1
    dt = 0.1
    D = 2
    #=== run the model =======================================================================
    d1 = Diffusion( X, T, dx, dt, D )
    d1.solve()
    d1.get_x_mean()
    d1.get_x2_mean()
    data = {'<x>': d1.x_means,
            '<x^2>':d1.x2_means}
    #=== plot the results ========================================================================
    import matplotlib.pyplot as plt
    #=== Plot <x>, <x^2> vs. t ========================================================================
    '''
    fig = plt.figure()
    for n, y in enumerate( data ):
        ax = fig.add_subplot( 2, 1, n + 1 )
        ax.scatter( range( d1.T ),
                    data[y],
                    c = 'r' )
    ax.set_xlim( 0, d1.T )
    #    ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
    plt.show()
    '''
    #=== Plot x positions with a few times diffusion========================================================================
    fig = plt.figure( figsize = ( 8.5, 11 ) )
    ax = fig.add_subplot( 1, 1, 1 )
#    for t in np.array( range( 20 ) ) ** 2:
    for t in np.linspace( 0, 25, num = 26 ) ** 2:
        ax.plot( range( d1.X ),
                 d1.u[:, t],
                 label = 't=%d' % t
                 )
    ax.set_ylim( ax.axis()[2] , ax.axis()[3] * 1.2 )
    ax.set_xlabel( 'X from 0 to % d' % d1.X )
    ax.set_ylabel( 'u( x, t )' )
    plt.legend()
    fig.suptitle( 'PHY260 Homework #6, Part 2: \nApproximation of the Diffusion Equation\ndt = %2.1f, dx = %2.1f, T from 0 to %d' % ( d1.dx, d1.dt, d1.T ) ,
                  ha = 'center', size = 'large' )
#    ax.set_xlim( 0, d1.X )
    # make axes symmetrical:
#    ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
    fig.savefig( 'figs/Diffusion_withoutcurves.pdf' )
