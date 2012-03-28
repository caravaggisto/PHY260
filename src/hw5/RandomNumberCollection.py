'''
Created on Mar 27, 2012

@author: res
'''

import numpy as np
import random
import scitools.std
from matplotlib.pyplot import plot, show
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from collections import defaultdict
from bisect import bisect_left

'''
a) Generate 1,000 random numbers, evenly distributed between 0 and 1.
Analyze the distribution of random numbers with increasing resolution,
by plotting the probability distribution of your random numbers, first
with 10 subdivisions, then 20, 50 and finally 100 subdivisions. Repeat the
exercise on a sample of 1,000,000 random numbers. [5 points]
'''
# sample size
N = 1000
y_values = np.asarray( range( N ) )
# next line makes us use the same random sample each time
# *** good for debugging!
random.seed( 1 )
random_sample = np.asarray( [random.uniform( 0, 1 ) for i in y_values] )

# make a plot figure
fig_samples = Figure( figsize = ( 8.5, 11 ) )
# Add title to top of figure
fig_samples.suptitle( '%d Samples of Uniformly Distributed ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
ax = fig_samples.add_subplot( 1, 1, 1 )
# plot the points
ax.scatter( random_sample, y_values, s = 2, color = 'blue' )
# add title to top of plot
ax.set_title( "%d Samples" % N, fontsize = 12 )
# set axis limits
ax.set_ylim( [0, N] )
ax.set_xlim( [0, 1] )
#save figure to a png file
canvas = FigureCanvas( fig_samples )
canvas.print_figure( 'samples_%d.png' % N, dpi = 500 )

fig_histograms = Figure( figsize = ( 8.5, 11 ) )
fig_histograms.suptitle( 'Histograms of %d Uniformly Distributed ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
for num_bins, subplot_id in zip( ( 10, 20, 100, 1000 ), ( 1, 2, 3, 4 ) ):
    ax = fig_histograms.add_subplot( 2, 2, subplot_id )
    # make a histogram of the sample
    bins = np.linspace( 0, 1 , num_bins )
    ax.hist( random_sample, bins, normed = 1, facecolor = 'green', alpha = 0.5, edgecolor = 'none' )
    # make nice labels, plot titles...
    ax.set_title( "%d Bins (width = %.3f)" % ( num_bins, random_sample[1] ), fontsize = 12 )
    ax.set_xlabel( "Value, between 0 and 1", fontsize = 10 )
    # add a nice grid.
    ax.grid( True, linestyle = '- ', color = '0.75' )
# save to a png file
canvas = FigureCanvas( fig_histograms )
canvas.print_figure( 'histograms_%d.png' % N, dpi = 500 )
