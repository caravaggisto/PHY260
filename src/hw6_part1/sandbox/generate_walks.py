'''
@author: willg

PHY260, Homework #6

This module generates data for problem 1, which writes:

1. 2D Random Walk
    Write a program to simulate a random walker in 2 dimensions, 
    taking steps of unit length in +/- x or +/- y direction on a 
    discrete square lattice.
'''
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

#===============================================================================
# Physical parameters for the problem 
#np.random.seed( 1 )    # set of randoms to use, for debugging.
Np = 1000                 # number of particles
Ns = 10000                # number of steps
#===============================================================================

#===============================================================================
# Logistical parameters
local_data_dir = 'data' # local directory in which npz files are stored.
local_figs_dir = 'figs' # and figures
data_path = get_data_directory( local_data_dir ) + '/' + "hw8"
figs_path = get_data_directory( local_figs_dir ) + '/' + "hw8"
#===============================================================================

if __name__ == '__main__':
    time_string = datetime.datetime.now().strftime( "%Y-%m-%d-%H:%M:%S" )
    # taken from http://stackoverflow.com/questions/6327498/
    walks = random_walk_2d( Np, Ns )
    data_path = data_path + '_' + time_string
    np.savez( data_path, walks = walks )
    # including x = x makes it easy to access the arrays later with file_handle['x'].
    # TODO: save metadata with this... like Ns, Np, time to completion, date+time stamp.
    # to do this, maybe PyTables is useful.


