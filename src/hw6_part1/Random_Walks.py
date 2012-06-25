'''
Created on Apr 22, 2012

@author: willg
'''
import numpy as np

# For debugging, to use same random set.
#np.random.seed( 1 )

class Random_Walks( object ):
    ''' Generates Np random walks of length Ns.
    Compute...
        <x> with get_x_means
        <x^2> with get_x2_means
    Methods to compute <x>, <x^2> are get_x_means, get
    '''
    def __init__( self, Np, Ns ):
        self.ns_values = np.arange( 0, Ns + 1 ) # 'Ns+1' because I want Ns *included* in the array.
        self.Np, self.Ns = Np, Ns
    def make_walks( self ):
        ''' computes and returns Np random walks with Ns steps each, organized thus:
            walks[p,s,i] = the ith coordinate *change* (-1,0, or 1) of the sth step for,j,i  the pth particle.
        e.g.
            walks[:,0:n,0] = the x-coord *changes* of all particles up to and including the nth step.
            walks[p,:,0].sum() = the final (after all steps) x-coord of the pth particle.
        '''
        Np, Ns = self.Np, self.Ns
        directions = {1:( 1, 0 ), 2:( 0, 1 ), 3:( -1, 0 ), 4:( 0, -1 )}
        random_ints = np.random.randint( 1, 5, Np * Ns )
        self.walks = np.zeros( ( Np * Ns, 2 ), dtype = int )
        for key, value in directions.iteritems():
            self.walks[random_ints == key] = value
        # resize is permanent; reshape is not.
        self.walks.resize( Np, Ns, 2 )
        return self.walks
    def get_x_means( self ):
        ''' computes and returns <x> '''
        ns_values, walks = self.ns_values, self.walks
        self.x_means = np.asarray( [walks[:, 0:n, 0].sum( axis = 1 ).mean() for n in ns_values] )
        return self.x_means
    def get_x2_means( self ):
        ''' computes and returns <x^2> '''
        ns_values, walks = self.ns_values, self.walks
        self.x2_means = np.asarray( [( walks[:, 0:n, 0].sum( axis = 1 ) ** 2 ).mean() for n in ns_values] )
        return self.x2_means
    def get_r2_mean( self ):
        ''' computes and returns <x^2> '''
        ns_values, walks = self.ns_values, self.walks
        self.r2_means = np.asarray( [( walks[:, 0:n, :].sum( axis = 1 ) ** 2 ).sum( axis = 1 ).mean() for n in ns_values] )
        return self.r2_means

if __name__ == '__main__':
    '''this code is for debugging the class, making sure 
    the methods are all working the way they're supposed to, etc. '''
    np.random.seed( 1 ) # use the same set of random numbers... forever!
    #===========================================================================
    Np = 2 # particles
    Ns = 4 # steps
    #===========================================================================
    w1 = Random_Walks( Np, Ns ) 
    w1.make_walks()
    w1.get_x_means()
    w1.get_x2_means()
    w1.get_r2_mean()
    data = {'<x>': w1.x_means,
            '<x^2>':w1.x2_means}
    # plot the results
    import matplotlib.pyplot as plt
    fig = plt.figure()
    for n, y in enumerate( data ):
        ax = fig.add_subplot( 2, 1, n + 1 )
        ax.scatter( w1.ns_values,
                    data[y],
                    c = 'r' )
        ax.set_xlim( 0, Ns )
    #    ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
#    plt.show()
