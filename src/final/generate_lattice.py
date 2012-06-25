'''
Created on May 1, 2012

@author: willg
'''
from Lattice import Lattice
from final import filepath_for_now, local_data_directory, \
    get_newest_file_by_type_in_dir, figs_directory
import pickle
import numpy as np

def make_stuff_for_part_1():
    ''' Find M vs. T '''
    #===============================================================================
    # Adjustable parameters for the problem 
    J = 1.5
    n = 50 #TODO: make this array-able... in B) we need to compare results for different n.
    n_steps = 20000
    kT_min = 1.0
    kT_max = 6.0
    kT_step = 0.1
    #===============================================================================
    kT_array = np.arange( kT_min, kT_max, kT_step )
    lattice = Lattice( J, n, n_steps, kT_array )
    lattice.run_metropolous_algorithm()
    
    fp = filepath_for_now( local_data_directory )
    print "Save to "
    print fp
    pickle.dump( lattice, open( fp + '.p' , "wb" ) )

def make_stuff_for_part_2():
    ''' Find C/N vs. n '''
    #===============================================================================
    # Adjustable parameters for the problem 
    J = 1.5
    n_list = [5, 10, 20, 30, 50, 75, 100, 200, 500] #TODO: make this array-able... in B) we need to compare results for different n.
    # TODO: comment out the following
#    n_list = [10]
    kT_min = 1.5
    kT_max = 4.5
    kT_step = 0.1
    #===============================================================================
    # make array of temps
    kT_array = np.arange( kT_min, kT_max + kT_step, kT_step )
    # make an empty array to store lattice simulations
    lattices = []
    for n in n_list:
        n_steps = int( ( n ** 2 ) * 100 )
        lattice = Lattice( J, n, n_steps, kT_array )
        lattice.run_metropolous_algorithm()
        lattices.append( lattice )
    # save the list of lattice simulations to a pickle:
    fp = filepath_for_now( local_data_directory )
    print "Save to "
    print fp
    pickle.dump( lattices, open( fp + '.p' , "wb" ) ) 

if __name__ == '__main__':
    make_stuff_for_part_2()

path_to_newest_pickle = get_newest_file_by_type_in_dir( '.p', 'data' )
path_for_new_fig_for_newest_pickle = figs_directory + '/' + path_to_newest_pickle.split( '/' )[-1][:-2]

