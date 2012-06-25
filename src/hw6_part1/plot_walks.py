'''
@author: willg

PHY260, Homework #6

This module plots results for problem 1, part a), which writes:

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
    title: "PHY260, Homework #6, 1.A: <x>, <x^2> with n.
    x-axis: "n"
    rasterized: False
'''
from hw6_part1.generate_walks import path_for_new_fig_for_newest_pickle, \
    path_to_newest_pickle
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.pyplot import legend

if __name__ == '__main__':
    # get the data from the Random_Walk object in a1_generate.py:
    random_walks = pickle.load( open( path_to_newest_pickle, "rb" ) )
    # make a fig for plotting, saving the data in `random_walks`
    fig = plt.figure( figsize = ( 8.5, 11 ) ) # 8.5" x 11"
    data = {'<X>' : random_walks.x_means,
            '<R^2>': random_walks.r2_means}
    ns_values = random_walks.ns_values
    colors = ['b', 'r']
    labels_xaxis = ['# of steps', '# of steps']
    labels_yaxis = ['Mean x-coordinate', 'Mean square of total distance from origin']
    symmetrical_yaxis = [True, False]
    add_regression = [False, True]
    for n, title in enumerate( data.keys() ):
        # Make a plot for both <x>, <x^2>
        ax = fig.add_subplot( 2, 1, n + 1 )
        ax.scatter( x = ns_values,
                    y = data[title],
                    s = 15, # size of dot
                    c = colors[n], # color
                    marker = 's'
                    )
        # Print title, axis labels, 
        ax.set_title( title )
        ax.set_xlabel( labels_xaxis[n] )
        ax.set_ylabel( labels_yaxis[n] )
        # Set axis limits. 
        ax.set_xlim( 0, random_walks.Ns )
        if symmetrical_yaxis[n]:
            # For the <x> plot, we want y=0 at the middle of the axis.
            ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
            # Looks nice and symmetrical that way.
        if add_regression[n]:
            # Regress, force y-intercept = 0.
            x = ns_values.reshape( ns_values.size, 1 )
            a, _, _, _ = np.linalg.lstsq( x, data[title] )
            # Calculate R-Squared
            r2 = ( np.corrcoef( data[title], a * ns_values )[0, 1] ) ** 2
            # Plot regression line, annotations with metadata. 
            ax.plot( x, a * x, 'y-' )
            x_annotation, y_annotation = 2.5, 102.5
            ax.text( x_annotation, y_annotation , 'Part b)' )
            y_annotation -= 10
            ax.text( x_annotation, y_annotation, 'Regression Coefficient = %8.6f' % a )
            y_annotation -= 10
            ax.text( x_annotation, y_annotation, 'R-squared = %8.6f' % r2 )
    # Give the whole figure a title.
    fig.suptitle( 'PHY260 Homework #6, Part 1: \n%d Random Walks with %d Steps' % ( random_walks.Np, random_walks.Ns ),
               ha = 'center', size = 'large' )
    legend()
    # Save the figure to a pdf file.
    fig.savefig( path_for_new_fig_for_newest_pickle + '.png' )


if False:
    ''' this code is from Random_Walks.py, 
    reproduced here as a reminder of some useful plotting tips,
    e.g. using a dict for managing subplots, etc '''
#    data = {'<x>': w1.x_means,
#            '<x^2>':w1.x2_means}
#    # plot the results
#    import matplotlib.pyplot as plt
#    fig = plt.figure()
#    for n, y in enumerate( data ):
#        ax = fig.add_subplot( 2, 1, n + 1 )
#        ax.scatter( w1.ns_values,
#                    data[y],
#                    c = 'r' )
#        ax.set_xlim( 0, Ns )
#    #    ax.set_ylim( -max( -ax.axis()[2], ax.axis()[3] ), max( -ax.axis()[2], ax.axis()[3] ) )
#    plt.show()
    pass
