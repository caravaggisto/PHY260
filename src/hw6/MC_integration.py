'''
Created on Apr 15, 2012

@author: willg
'''

data_filename_root = "MCint_data"
data_dir = 'data' # local directory in which npz files are stored.
figs_dir = 'figs' # and figures

import numpy as np
import os
import datetime

def montecarlo_integration_v1( f, a, b, n ):
    ''' monte carlo integration of
        f = callable function
        a = start coordinate (an int)
        b = end coordinate
        n = number of iterations for the approximation (resolution of the MC integral)
    returns I = integral approximation.
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
    returns k_vals, I_vals
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

def f1( x ):
    return 2 + 3 * x

if __name__ == '__main__':
    time_string = datetime.datetime.now().strftime( "%Y-%m-%d-%H:%M:%S" )
    # taken from http://stackoverflow.com/questions/6327498/
    a, b = 1, 2
    n = int( 1e7 )
    k, I = montecarlo_integration_v3( f1, a, b, n, N = int( n / 1000 ) )
    data_filepath = get_data_directory( data_dir ) + '/' + data_filename_root + '_' + time_string
    np.savez( data_filepath, k = k, I = I )
    # including x = x makes it easy to access the arrays later with file_handle['x'].
