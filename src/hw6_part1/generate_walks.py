'''
@author: willg

PHY260, Homework #6

This module produces results for problem 1, part a), which writes:

1. a) Plot <x_n> and <(x_n)^2> up to n = 100 by averaging over at 
    least 1e4 different walks for each n > 3. 

--> Ns = n elem of {0,1,2,...,100}
--> Np = 1e4 = 10000
in this problem (1.a) we want two SCATTER plots on the same figure:
    1. <x> vs. n
        a. y-axis: "<x>"
        a. y-lims: plus or minus 1.2 times the largest magnitude data point.
    2. <x^2> vs. n
        a. y-axis: "<x^2>"
other notes:
    name: "PHY260, Homework #6, 1.A: <x>, <x^2> with n.
    x-axis: "n"
    rasterized: False
'''

from file_manager import local_data_directory, get_newest_file_by_type_in_dir, \
    filepath_for_now
import numpy as np
import pickle

class Random_Walks( object ):
    ''' Generates Np random walks of length Ns.
    Compute...
        <x> with get_x_mean
        <x^2> with get_x2_mean
    Methods to compute <x>, <x^2> are get_x_mean, get
    '''
    def __init__( self, Np, Ns ):
        self.ns_values = np.arange( 1, Ns + 1 ) # 'Ns+1' because I want Ns *included* in the array.
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
    def get_x_mean( self ):
        ''' computes and returns <x> '''
        n_values, walks = self.ns_values, self.walks
        self.x_means = np.asarray( [walks[:, 0:n, 0].mean() for n in n_values] )
        return self.x_means
    def get_x2_mean( self ):
        ''' computes and returns <x^2> '''
        n_values, walks = self.ns_values, self.walks
        self.x2_means = np.asarray( [np.square( walks[:, 0:n, 0] ).mean() for n in n_values] )
        return self.x2_means

if __name__ == '__main__':
    #===============================================================================
    # Adjustable parameters for the problem 
    #TODO: make these 10000 and 100, respectively.
    Np = 100
    Ns = 10
    #===============================================================================
    random_walks = Random_Walks( Np, Ns )
    random_walks.make_walks()
    random_walks.get_x_mean()
    random_walks.get_x2_mean()
    pickle.dump( random_walks, open( filepath_for_now( local_data_directory ) + '.p' , "wb" ) )

newest_pickle_file = get_newest_file_by_type_in_dir( '.p' , local_data_directory )
newest_data_object = pickle.load( open( newest_pickle_file, "rb" ) )
#data_path = data_filename_root + '_' + time_string + '.p'
