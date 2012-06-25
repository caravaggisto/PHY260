'''
Created on Apr 15, 2012

@author: willg

WARNING: don't run before actually generating data with MC_integration.
'''

from MC_integration import data_dir, figs_dir
import numpy as np
import matplotlib.pyplot as plt
import os
from Image import RASTERIZE

def get_newest_file_in_data_dir():
    # taken from http://ubuntuforums.org/showthread.php?t=1526010
    os.chdir( data_dir )
    # change into the data directory
    npz_filelist = [npz_file for npz_file in os.listdir( os.getcwd() ) if npz_file.endswith( '.npz' ) ]
    # make a list of all files and folders in the current working directory
    npz_filelist = filter( lambda x: not os.path.isdir( x ), npz_filelist )
    # get rid of folders in the list. we only want files.
    newest = max( npz_filelist, key = lambda x: os.stat( x ).st_mtime )
    # get newest file in the list
    os.chdir( '..' )
    # change back to the parent directory.
    return newest

def get_data( np_file_path ):
    data = np.load( np_file_path )
    k, I = data['k'], data['I']
    error = 6.5 - I # 6.5 = analytical solution, I is np array of approximations
    return k, error

def make_data_figure( x, y ):
    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.plot( x, y )
    return fig

if __name__ == '__main__':
    npzfilename = get_newest_file_in_data_dir()
    k, error = get_data( os.getcwd() + '/' + data_dir + '/' + npzfilename )
    # plot the results for the solution as a function of n:
    fig = make_data_figure( k, error )
    fig.savefig( os.getcwd() + '/' + figs_dir + '/' + npzfilename[:-3] + 'pdf' )
