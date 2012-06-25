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

from Random_Walks import Random_Walks
from hw6_part1 import filepath_for_now, local_data_directory, \
    get_newest_file_by_type_in_dir, figs_directory
import pickle

def make_walks():
    #===============================================================================
    # Adjustable parameters for the problem 
    Np = 10000
    Ns = 100
    #===============================================================================
    random_walks = Random_Walks( Np, Ns )
    random_walks.make_walks()
    random_walks.get_x_means()
    random_walks.get_x2_means()
    random_walks.get_r2_mean()
    fp = filepath_for_now( local_data_directory )
    print "Save to "
    print fp
    pickle.dump( random_walks, open( fp + '.p' , "wb" ) )

if __name__ == '__main__':
    make_walks()

path_to_newest_pickle = get_newest_file_by_type_in_dir( '.p', 'data' )
path_for_new_fig_for_newest_pickle = figs_directory + '/' + path_to_newest_pickle.split( '/' )[-1][:-2]

