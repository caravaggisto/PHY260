'''
Created on Feb 3, 2012

@author: res
'''
from PHY260_Grisaitis_homework2_Radioactive_Decay import Radioactive_Decay
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants.constants import N_A
from matplotlib.ticker import ScalarFormatter

class Data_Plot_Figure( object ):
    '''
    classdocs
    '''
    
    def __init__( self, figure_number = 1, \
                  figure_title = 'New Figure', \
                  x_axis_label = 'x-axis', \
                  y_axis_label = 'y-axis', \
#                  figure_title = "Approximations and Exact Solution for Radioactive Decay of Carbon-14", \
#                  x_axis_label = "$t$ in years", \
#                  y_axis_label = "$N(t)$, number of atoms remaining", \
                  figure_size = ( 11.0, 8.5 ), \
                  axes_dimensions = [0.05, 0.1, 0.9, 0.8], \
                  legend_position = 'lower left', \
                  background_color = 'white', \
                  ):
        '''
        Constructor
        '''
#        TODO: Do I even need a figure number if I'm enclosing all of this stuff in an object?
        self.figure_number = figure_number
#        self.figure_size = figure_size
#        self.axes_dimensions = axes_dimensions
        self.time_arrays = []
        self.data_arrays = []
        self.legend_labels = []
        self.plots = []
        # make a new figure, callable with self.plot. Specify size.
        self.plot = plt.figure( self.figure_number, figsize = figure_size, )
        #TODO: Does the following line work? 
        self.axes = self.plot.add_axes( axes_dimensions )
#        self.axes = plt.figure( self.figure_number ).add_axes( axes_dimensions )
        self.axes.grid( True )        # include a grid
        # Set titles, labels, white background
        self.axes.set_title( figure_title )
        self.axes.set_xlabel( x_axis_label )
        self.axes.set_ylabel( y_axis_label )
        #TODO: Does the following line work? The two following do, but it's verbose and weird.
        self.plot.patch.set_facecolor( background_color )
#        self.rect = plt.figure( self.figure_number ).patch
#        self.rect.set_facecolor( background_color )
        
        #TODO: Do I need these variables? (self.figure_title...)
#        self.figure_title = figure_title
#        self.x_axis_label = x_axis_label
#        self.y_axis_label = y_axis_label
        
        # Make axis tick numbers look nice
        # TODO: do I have to do this anywhere else, say in update_plot_axis()?
        niceMathTextForm = ScalarFormatter( useMathText = True )
#        self.axes.xaxis.set_major_formatter( niceMathTextForm )
        self.axes.yaxis.set_major_formatter( niceMathTextForm )
    
    def update_figure_data( self, time_array, data_array, data_name ):
        self.time_arrays.append( time_array )
        self.data_arrays.append( data_array )
        self.legend_labels.append( data_name )
    
    def update_plot_legend_and_axes( self ):
        self.plot.legend( self.plots, self.legend_labels, self.legend_position )
#        syntax: self.plot.ylim( lower bound, upper bound )
        self.plt.ylim( min( 0, min( self.data_arrays ) ), max( self.data_arrays ) * 1.2 * N_A )
        self.plt.xlim( min( self.time_arrays ), max( self.time_arrays ) )
    
    def plot_data_scatter( self, time_array, data_array, data_name, marker_color, marker_type, marker_size, marker_border_size ):
        self.update_figure_data( time_array, data_array, data_name )
        self.plots.append( self.axes.scatter( time_array, \
#                                    mol_to_particles( data_array ), \
                                    data_array, \
                                    c = marker_color, \
                                    marker = marker_type, \
                                    s = marker_size, \
                                    lw = marker_border_size \
                                    ) )
        self.update_plot_legend_and_axes()
        return
    
    def plot_data_line( self, time_array, data_array, data_name, line_color, line_type, line_width ):
        self.update_figure_data( time_array, data_array, data_name )
        self.plots.append( self.axes.plot( time_array, \
#                                    mol_to_particles( data_array ), \
                                    data_array, \
                                    line_color + line_type, \
                                    line_width \
                                    ) )        
        self.update_plot_legend_and_axes()
    
    def save_figure_to_file( self, file_directory, file_name_base, file_types ):
        ''' save figure to file of each type in file_types. '''
        for file_type in file_types:
            plt.savefig( file_directory + file_name_base + "." + file_type, facecolor = self.plot.get_facecolor(), edgecolor = 'none' )
    
    def show_figure( self ):
        ''' show window with figure '''
        self.plt.show()
    
    def set_figure_title( self, figure_title ):
        # TODO: Do I need this method?
        self.axes.set_title( figure_title )
    
    def set_x_axis_label( self, x_axis_label ):
        # TODO: Do I need this method?
        self.axes.set_xlabel( x_axis_label )
    
    def set_y_axis_label( self, y_axis_label ):
        # TODO: Do I need this method?
        self.axes.set_ylabel( y_axis_label )
    
