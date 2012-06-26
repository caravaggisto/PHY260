'''
Created on Jan 30, 2012

@author: res
'''

from PHY260_Grisaitis_homework2_Radioactive_Decay import Radioactive_Decay
from matplotlib.ticker import ScalarFormatter
from scipy.constants.constants import N_A
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # file system variables
    file_directory = "/tmp/"
#    TODO: edit the following line, implement differently when processing data for different parts of the assignment.
    file_name_base = "PHY260_Grisaitis_homework2_" + "PartB"
    figure_file_extensions = ["pdf", "eps", "png"]
    
    # physical contants
    c14_molar_mass = 14.003241  # grams per mole (source: http://en.wikipedia.org/wiki/Carbon-14 )
    c14_half_life = 5700.0      # years
    
    # modeling parameters
    # time
    start = 0
    stop = 20000
    partB_time_step_widths = [10, 100]
    # initial condition
    initial_mass = 10.0 ** -14    # kg
    
    # helper function    
    def mol_to_particles( mol ):
        return mol * N_A
    
    '''Part B: '''
    # initialize a radioactive decay model :)
    decay_simulation = Radioactive_Decay( initial_mass, c14_molar_mass, c14_half_life, )
    initial_moles = decay_simulation.N0
    # make lists for results and time arrays, for graphing later
    results = []
    time_arrays = []
    # simulation parameters
    time_step_widths = [10, 100]
    # run Euler approximations for time widths 10 and 100 years.
    for k in range( 2 ):
        # add 1 for including the end point
        num = 1 + ( stop - start ) / time_step_widths[k]
        # include endpoint T. E.g.: [0,1,2,...,T]
        time_arrays.append( np.linspace( start, stop, num, endpoint = True ) )
        t, N = decay_simulation.approximation_Forward_Euler( time_arrays[k] )
        #TODO: during cleanup, remove this plotting code:
#            plt.figure( k, figsize = ( 11.0, 8.5 ) )
#            axes = plt.figure( k ).add_axes( [0.1, 0.2, 0.8, 0.7] )
#            axes.plot( t, N, 'b-', 1 )
#            plt.ylim( 0., initial_moles * 1.2 )
#            for file_type in ["pdf", "eps", "png"]:
#                plt.savefig( file_directory + str( k ) + "." + file_type )
#            plt.show()
        results.append( N )
    # store N(t)'s for the analytical (exact) solution 
    t, N_exact = decay_simulation.analytical_solution_array( time_arrays[0] )
    results.append( N_exact )
    time_arrays.append( t )
    
    
    '''plot stuff with matplotlib.pyplot'''
    # make a new figure, callable with plt.figure(2). Specify size.
    plt.figure( 2, figsize = ( 11.0, 8.5 ) )
    axes = plt.figure( 2 ).add_axes( [0.05, 0.1, 0.9, 0.8] )
    # include a grid
    axes.grid( True )
    point_labels = ['ko', 'k+', 'ro']
    point_sizes = [30, 100, 4]
    num_of_graphs = len( point_labels )
    line_widths = [0, 1, 0]
    #Next 3 lines set titles and labels and whitebackground
    axes.set_title( "Approximations and Exact Solution for Radioactive Decay of Carbon-14" )
    axes.set_xlabel( "$t$ in years" )
    axes.set_ylabel( "$N(t)$, number of atoms remaining" )
    rect = plt.figure( 2 ).patch
    rect.set_facecolor( 'white' )
    # Next three lines make the axis tick numbers look nice
    niceMathTextForm = ScalarFormatter( useMathText = True )
#        axes.xaxis.set_major_formatter( niceMathTextForm )
    axes.yaxis.set_major_formatter( niceMathTextForm )
    # make plot objects
    plots = []
    for k in range( num_of_graphs ):
#            plots.append( axes.plot( time_arrays[k], mol_to_particles( results[k] ), point_labels[k], line_widths[k] ) )
        plots.append( axes.scatter( time_arrays[k], \
                                    mol_to_particles( results[k] ), \
                                    c = point_labels[k][0], \
                                    marker = point_labels[k][1], \
                                    s = point_sizes[k], \
                                    lw = line_widths[k] \
                                    ) )

#            plots.append( axes.plot( time_arrays[k], results[k], point_labels[k], line_widths[k] ) )
#        plt.legend( plots, ( "$\\sin(x)$", "$\\sin(x)/x$" ), 'lower center' )
    legend_labels = ["Approximation; Steps of 10 years", "Approximation; Steps of 100 years", "Exact solution"]
    legend_position = 'lower left'
    plt.legend( plots, legend_labels, legend_position )
#        plt.ylim( 0., initial_moles * 1.2 )
    plt.ylim( 0., initial_moles * 1.2 * N_A )
    plt.xlim( start, stop )
    # save files
    for file_type in figure_file_extensions:
        plt.savefig( file_directory + file_name_base + "part_B_plots" + "." + file_type, facecolor = plt.figure( 2 ).get_facecolor(), edgecolor = 'none' )
#        TODO: uncomment the following statement.
    plt.show()
    
    ''' Part C '''
    
    


'''
#####################################
Ways in which this could be improved:
#####################################
1. Move the graphing and figure management into another class, 
    - probably in the ODE solver class. There, I could manage what 
    figure were being manipulated, return the figure objects back, etc.
2. Move the question parts into classes.
    - Make class "Part". Each part would have some simulation and 
    graphing in common, but different parameters to try out, etc.

'''
