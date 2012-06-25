'''
@author: willg

PHY260, Homework #6

This module plots results for problem 1, part a), which writes:

1. a) Plot <x_n> and <(x_n)^2> up to n = 100 by averaging over at 
    least 1e4 different walks for each n > 3. 
'''

from generate_walks import local_figs_dir, Ns, Np
from retrieve_walks import walks, npzfilename
import numpy as np
import matplotlib.pyplot as plt
import os

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
    # make and save the figure...
    fig = make_data_figure( walks[:, 0], walks[:, 1] )
    fig.savefig( os.getcwd() + '/' + local_figs_dir + '/' + npzfilename[:-3] + 'pdf' )




