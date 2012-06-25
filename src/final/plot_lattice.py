'''
@author: willg

PHY260, Final Project

This module plots the results a Metropolous algorithm for a lattice
of n x n particles. It's a type of Ising model.

a) Plot M vs. T
--> get array of temperature values
--> load the data generated for those temps.
Determine the critical temperature T_C
--> eyeball the value and annotate.

b) Calculate the specific heat per spin (C/N) for the following lattice sizes:
    n in [5,10,20,30,40,50,75,100,200,500]
using the fluctuation-dissipation theorem:
    C = (1/k_B)*(dE/T)^2
and verify that
    Cmax / N ~ log(n)
--> make array of log(n)
--> determine Cmax for each n

    title: "PHY260, Final Exam, A: Total Magnetization as a function of k_B*T.
    x-axis: "k_B*T"
    y-axis: "Total Magnetization = M = N<s>"
    rasterized: False
'''
from final.generate_lattice import path_for_new_fig_for_newest_pickle, \
    path_to_newest_pickle
import matplotlib.pyplot as plt
import numpy as np
import pickle

def plot_stuff_for_part1():
    # get the data from the Lattice object:
    lattice = pickle.load( open( path_to_newest_pickle, "rb" ) )
    # make a fig for plotting, saving the data in `lattice`
    fig = plt.figure( figsize = ( 7, 5 ) )
    data = {'M vs kT' : lattice.M }
    labels_xaxis = ['temperature (kT)', '']
    labels_yaxis = ['Total Mag', '']
    symmetrical_yaxis = [True, False]
    for n, y in enumerate( data ):
        ax = fig.add_subplot( len( data ), 1, n + 1 )
        ax.scatter( lattice.kT_array,
                    data[y],
                    c = 'r' )
        ax.set_xlabel( labels_xaxis[n] )
        ax.set_ylabel( labels_yaxis[n] )
        ax.set_xlim( lattice.kT_array.min(), lattice.kT_array.max() )
        if symmetrical_yaxis[n]:
            # y=0 at the middle of the axis.
            ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
            # Looks nice and symmetrical now.
    # Give the whole figure a title.
    fig.suptitle( 'PHY260 Final Exam, Part 1: \n%d x %d Lattice with %d Relaxations' % ( lattice.n, lattice.n, lattice.n_steps ),
               ha = 'center', size = 'large' )
    # Save the figure to a pdf file.
    fig.savefig( path_for_new_fig_for_newest_pickle + 'M' + '.png' )

def part2_all_specific_heats( lattices ):
    # Parameters for the plot figure:
    labels_xaxis = 'temperature (kT)'
    labels_yaxis = 'Specific Heat (per lattice point)'
    # make a fig for plotting, saving the data in `lattice`
    fig = plt.figure( figsize = ( 7, 5 ) )
#=== Subplot for Specific Heat============================================================================
    ax = fig.add_subplot( 2, 1, 1 )
    for n in range( len( lattices ) ):
        print n
        lattice = lattices[n]
        ax.plot( lattice.kT_array,
                 lattice.C,
                 label = ( 'n = %d,\nrelaxes = %d' % ( lattice.n, lattice.n_steps ) )
                 )
    # make legend, set axis limits
    plt.legend()
    ax.set_xlabel( labels_xaxis )
    ax.set_ylabel( labels_yaxis )
    ax.set_xlim( lattice.kT_array.min(), lattice.kT_array.max() )
#=== Subplot for Magnetization============================================================================
    ax = fig.add_subplot( 2, 1, 2 )
    for n in range( len( lattices ) ):
        print n
        lattice = lattices[n]
        ax.plot( lattice.kT_array,
                 lattice.M,
                 label = ( 'n = %d,\nrelaxes = %d' % ( lattice.n, lattice.n_steps ) )
                 )
    # make legend, set axis limits
    plt.legend()
    ax.set_xlabel( labels_xaxis )
    ax.set_ylabel( 'M' )
    ax.set_xlim( lattice.kT_array.min(), lattice.kT_array.max() )
    # Give the whole figure a title.
    fig.suptitle( 'PHY260 Final Exam, Part 2: \nSpecific Heat and Magnetization for n = %d' % lattice.n,
               ha = 'center', size = 'large' )
    # Save the figure to a pdf file.
    fig.savefig( path_for_new_fig_for_newest_pickle + '_C' + '.png' )

def part2_plot_C_and_M_together( lattices ):
    # Parameters for the plot figure:
    labels_xaxis = 'temperature (kT)'
    # make a fig for plotting, saving the data in `lattice`
    fig = plt.figure( figsize = ( 7, 5 ) )
#=== Subplot for Specific Heat============================================================================
    ax = fig.add_subplot( 2, 1, 1 )
    for n in range( len( lattices ) ):
        print n
        lattice = lattices[n]
        ax.plot( lattice.kT_array,
                 lattice.C,
                 label = ( 'n = %d,\nrelaxes = %d' % ( lattice.n, lattice.n_steps ) )
                 )
    # make legend, set axis limits
    plt.legend()
    ax.set_xlabel( labels_xaxis )
    ax.set_ylabel( 'C/N' )
    ax.set_xlim( lattice.kT_array.min(), lattice.kT_array.max() )
#=== Subplot for Magnetization============================================================================
    ax = fig.add_subplot( 2, 1, 2 )
    for n in range( len( lattices ) ):
        print n
        lattice = lattices[n]
        ax.plot( lattice.kT_array,
                 lattice.M,
                 label = ( 'n = %d,\nrelaxes = %d' % ( lattice.n, lattice.n_steps ) )
                 )
    # make legend, set axis limits
    plt.legend()
    ax.set_xlabel( labels_xaxis )
    ax.set_ylabel( 'M' )
    ax.set_xlim( lattice.kT_array.min(), lattice.kT_array.max() )
    # Give the whole figure a title.
    fig.suptitle( 'PHY260 Final Exam, Part 2: \nSpecific Heat and Magnetization for n = %d' % lattice.n,
               ha = 'center', size = 'large' )
    # Save the figure to a pdf file.
    fig.savefig( path_for_new_fig_for_newest_pickle + '_C&M' + '.png' )

if __name__ == '__main__':
    # get the data from the list of Lattice objects:
    print "Loading lattices in", path_to_newest_pickle, "..."
    lattices = pickle.load( open( path_to_newest_pickle, "rb" ) )
    print "Lattices loaded"
    part2_plot_C_and_M_together( lattices )
    
