'''
@author: willg
'''

data_filename_root = "randwalk_data"
data_dir = 'data' # local directory in which npz files are stored.
figs_dir = 'figs' # and figures

import datetime
import numpy as np
import os

def get_data_directory( local_dir ):
    ''' returns directory in which data is to be stored.
    if the directory doesn't exist, it creates it. 
    code taken from http://stackoverflow.com/questions/273192/
    '''
    path_to_data_dir = os.getcwd() + '/' + local_dir
    if not os.path.exists( path_to_data_dir ):
        os.makedirs( path_to_data_dir )
    # code taken from http://stackoverflow.com/questions/273192/
    return path_to_data_dir

'''
1. 2D Random Walk
write a program to simulate a random walker in 2 dimensions, 
taking steps of unit length in +/- x or +/- y direction on a 
discrete square lattice.
'''

#===============================================================================
# Parameters for the problem 
#np.random.seed( 1 )    # set of randoms to use, for debugging.
Np = 1000                 # number of particles
Ns = 10000                # number of steps
#===============================================================================

def random_walk_2d( Np, Ns ):
    ''' Generates data for Part 1.
    Np = # particles
    Ns = # steps
    '''
    positions = np.zeros( ( Np, 2 ) )
    directions = {1:( 1, 0 ), 2:( 0, 1 ), 3:( -1, 0 ), 4:( 0, -1 )}
    # a dictionary for handling grid movement.
    for step in range( Ns ):
        for particle in range( Np ):
            random_direction = np.random.randint( 1, 5 )
            # generates random int in [1,5).
            positions[particle] += directions[random_direction]
            # progress along the random step
            # alternative implementation: 
            #random_direction = random.choice(directions)
            #positions[particle] += random_direction
    return positions

if __name__ == '__main__':
    time_string = datetime.datetime.now().strftime( "%Y-%m-%d-%H:%M:%S" )
    # taken from http://stackoverflow.com/questions/6327498/
    positions = random_walk_2d( Np, Ns )
    data_filepath = get_data_directory( data_dir ) + '/' + data_filename_root + '_' + time_string
    np.savez( data_filepath, positions = positions )
    # including x = x makes it easy to access the arrays later with file_handle['x'].


#positions = random_walk_2d( Np, Ns )
#print positions[:, 0], positions[:, 1]
#print 2 * positions.std( axis = 0 ) ** 2 # approximately Ns!


'''

a)
plot <x_n> and <(x_n)^2> up to n = 100 by averaging over at 
least 1e4 different walks for each n > 3. 

b)
show that the motion is diffusive, i.e. that the mean square 
distance from the starting point r^2 is proportional to t 
(with t being the time, i.e. t ~ n) and determine the value 
of the diffusion constant (a simple "eyeball" fit to your 
numerical data is sufficient).
'''



