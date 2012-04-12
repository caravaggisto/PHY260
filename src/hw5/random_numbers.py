'''
Created on Mar 27, 2012

@author: William C Grisaitis

'''

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib import mlab
from matplotlib.backends.backend_pdf import PdfPages

'''
a) Generate 1,000 random numbers, evenly distributed between 0 and 1.
Analyze the distribution of random numbers with increasing resolution,
by plotting the probability distribution of your random numbers, first
with 10 subdivisions, then 20, 50 and finally 100 subdivisions. Repeat the
exercise on a sample of 1,000,000 random numbers. [5 points]
'''

def generate_nice_uniform( pdf_handler, N ):
    y_values = np.arange( N )
    # next line makes us use the same random sample each time
    # *** good for debugging!
    np.random.seed( 1 )
#    random_sample = np.asarray( [random.uniform( 0, 1 ) for i in y_values] )
    random_sample = np.random.uniform( 0, 1, N )
    
    # make a plot figure
    fig_samples = Figure( figsize = ( 8.5, 11 ) )
    # Add title to top of figure
    fig_samples.suptitle( 'Part A: %d Samples of Uniform ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
    ax = fig_samples.add_subplot( 1, 1, 1 )
    # plot the points
    ax.scatter( random_sample, y_values, s = 1, color = 'blue', lw = 0, rasterized = True )
    # add title to top of plot
    ax.set_title( "%d Samples" % N, fontsize = 12 )
    # set axis limits
    ax.set_ylim( [0, N] )
    ax.set_xlim( [0, 1] )
    #save figure to a png file
    canvas = FigureCanvas( fig_samples )
    canvas.print_figure( 'a_samples_%d.png' % N, dpi = 500 )
    # save figure to a pdf file
    fig_samples.savefig( pdf_handler, format = 'pdf' )
    
    fig_histograms = Figure( figsize = ( 8.5, 11 ) )
    fig_histograms.suptitle( 'Part A: Histograms of %d Uniform ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
    for num_bins, subplot_id in zip( ( 10, 20, 100, 1000 ), ( 1, 2, 3, 4 ) ):
        ax = fig_histograms.add_subplot( 2, 2, subplot_id )
        # make a histogram of the sample
        bins = np.linspace( 0, 1 , num_bins )
        ax.hist( random_sample, bins, normed = 1, facecolor = 'green', alpha = 0.5, edgecolor = 'none' )
        # make nice labels, plot titles...
        ax.set_title( "%d Bins (width = %.3f)" % ( num_bins, bins[1] ), fontsize = 12 )
        ax.set_xlabel( "Value, between 0 and 1", fontsize = 10 )
        # add a nice grid.
        ax.grid( True, linestyle = '- ', color = '0.75' )
    # save to a png file
    canvas = FigureCanvas( fig_histograms )
    canvas.print_figure( 'a_histograms_%d.png' % N, dpi = 500 )
    fig_histograms.savefig( pdf_handler, format = 'pdf' )

#######################################################################################

"""
b) Do the same as in part (a), but now with a sample of random numbers 
distributed according to a Gaussian distributuon with sigma = 1.0.
Verify your algorithm for generating random numbers by overlaying a 
Gaussian onto the plots of the distributions you generate. [5 points]
"""
def generate_nice_normal( pdf_handler, mu, sigma, N ):
    y_values = np.arange( N )
    # next line makes us use the same random sample each time
    # *** good for debugging!
    np.random.seed( 1 )
#    random_sample = np.asarray( [random.normalvariate( mu, sigma ) for i in y_values] )
    random_sample = np.random.normal( mu, sigma, N )
    # make a plot figure
    fig_samples = Figure( figsize = ( 8.5, 11 ) )
    # Add title to top of figure
    fig_samples.suptitle( 'Part B: %d Samples of Normal ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
    ax = fig_samples.add_subplot( 1, 1, 1 )
    # plot the points
    ax.scatter( random_sample, y_values, s = 1, color = 'blue', lw = 0, rasterized = True )
    # add title to top of plot
    ax.set_title( "%d Samples" % N, fontsize = 12 )
    # set axis limits
    ax.set_ylim( [0, N] )
    #ax.set_xlim( [0, 1] )
    #save figure to a png file
    canvas = FigureCanvas( fig_samples )
    canvas.print_figure( 'b_samples_%d.png' % N, dpi = 500 )
    fig_samples.savefig( pdf_handler, format = 'pdf' )
    
    fig_histograms = Figure( figsize = ( 8.5, 11 ) )
    fig_histograms.suptitle( 'Part B: Histograms of %d Normal ~ [0,1] Random Variables' % ( N, ), fontsize = 14 )
    x_vals = np.linspace( random_sample.min(), random_sample.max(), 1000 )
    norm_pdf = mlab.normpdf( x_vals, mu, sigma )
    for num_bins, subplot_id in zip( ( 10, 20, 100, 1000 ), ( 1, 2, 3, 4 ) ):
        ax = fig_histograms.add_subplot( 2, 2, subplot_id )
        # make a histogram of the sample
        ax.hist( random_sample, bins = num_bins, normed = 1, facecolor = 'green', alpha = 0.5, edgecolor = 'none' )
        # overlay a normal probability density function
        ax.plot( x_vals, norm_pdf, 'r--', linewidth = 1 )
        # make nice labels, plot titles...
        ax.set_title( "%d Bins" % ( num_bins ), fontsize = 12 )
        ax.set_xlabel( "Random Value", fontsize = 10 )
        # add a nice grid.
        ax.grid( True, linestyle = '- ', color = '0.75' )
    # save to a png file
    canvas = FigureCanvas( fig_histograms )
    canvas.print_figure( 'b_histograms_%d.png' % N, dpi = 500 )
    fig_histograms.savefig( pdf_handler, format = 'pdf' )

#######################################################################################
# Generate files:

if __name__ == '__main__':
    pdf_thing = PdfPages( 'PHY260_WilliamGrisaitis_homework#5.pdf' )
    for N in ( 1000, 1000000 ):
        generate_nice_uniform( pdf_thing, N )
        mu, sigma = 0, 1
        generate_nice_normal( pdf_thing, mu, sigma, N )
    # close the pdf file with all the figures.
    pdf_thing.close()
