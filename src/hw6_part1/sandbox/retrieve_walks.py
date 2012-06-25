'''
Created on Apr 15, 2012

@author: willg

WARNING: don't run before actually generating data with random_walk_data.
'''

from random_walk_data import data_dir, figs_dir, Ns, Np
import numpy as np
import matplotlib.pyplot as plt
import os

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
    positions = data['positions']
    return positions

def make_data_figure( x, y ):
    fig = plt.figure()
    ax = fig.add_subplot( 111 )
    ax.plot( x, y, 'ko', rasterized = True )
    ax.set_title( 'PHY260 Homework #6, Part 1: Plot of %d particles after %d steps' % ( Np, Ns ) )
    lim = max( abs( np.asarray( plt.axis() ) ) ) * 1.2
    plt.axis( [-lim, lim, -lim, lim] )
    # make sure the plot axes include all points, with a 20% margin from the closest point.
    plt.show()
    return fig

if __name__ == '__main__':
    print "hi"
    npzfilename = get_newest_file_in_data_dir()
    print "hi"
    positions = get_data( os.getcwd() + '/' + data_dir + '/' + npzfilename )
    print "hi"
    # make and save the figure...
    fig = make_data_figure( positions[:, 0], positions[:, 1] )
    print "hi"
    fig.savefig( os.getcwd() + '/' + figs_dir + '/' + npzfilename[:-3] + 'pdf' )
