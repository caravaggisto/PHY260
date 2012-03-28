'''
Created on Mar 1, 2012

@author: res
'''

from matplotlib.pyplot import figtext, figure, plot, sca, axes, axis
from numpy import asarray

class Grapher( object ):
    '''
    This class manages a figure with a single graph of 
    multiple line graphs. Initialize it with:
        title, subtitle, x_label, y_label, legend_loc
    Before you save any data, be sure to run these methods first:
        add_data( t, u, name )    array, array, string
    '''
    def __init__( self, title, subtitle, x_label, y_label, legend_loc ):
        '''
        Constructor
        '''
        self.title, self.subtitle, self.x_label, self.y_label, self.legend_loc = title, subtitle, x_label, y_label, legend_loc
        self.data_labels = []
        self.colors = ['r', 'y', 'g', 'b', 'c', 'm', 'k']
        self.num_colors = len( self.colors )
        self.k = 0      # for cycling through colors
        self.text_size = 10
        self.line_width = 1.0
        # Turn on the grid on the plot
        self.fig = figure( figsize = ( 7.5, 7.5 ) )
        self.ax = self.fig.add_subplot( 111 )
        self.ax.grid( True )
        self.fig.patch.set_facecolor( 'white' )
#        title( 'Simulated Trajectories of a Golf Ball (%s Method)' % method )
        self.ax.set_title( title )
        self.ax.text( .5, .85, subtitle, fontsize = self.text_size, ha = 'center' )
        self.ax.set_xlabel( x_label )
        self.ax.set_ylabel( y_label )
    
    def add_data( self, time_points, u_points, name ):
        ''' initialize with (array, array, string) '''
        #TODO: see if this works with label = name
        try:
            self.ax.plot( time_points, u_points, color = self.colors[self.k % self.num_colors], label = name, lw = self.line_width )
        except ValueError:
            self.ax.plot( time_points, u_points, color = self.colors[self.k % self.num_colors], lw = self.line_width )
#        self.data_labels.append( name )
        #TODO: make way of updating the legend. Does this work?
        self.k += 1 # update the color counter
        self.set_legend()
    
    def add_data_scatter( self, x, y, name ):
        ''' initialize with (array, array, string) '''
        self.ax.scatter( x, y, s = 5, c = 'r', marker = 'o', label = name, lw = 0 )
        self.set_legend()
    
    def set_legend( self ):
        """ Because I keep forgetting,
        String     Number
        upper right     1
        upper left     2
        lower left     3
        lower right     4
        right     5
        center left     6
        center right     7
        lower center     8
        upper center     9
        center     10
        """
        self.ax.legend( loc = self.legend_loc )
#        self.ax.relim()    # doesn't seem to help
        self.scale_axis_limits( 1, 1, 1, 1.3 )
    
    def set_axis_limits( self, xmin = None, xmax = None, ymin = None, ymax = None ):
        sca( self.ax )
        [xmin_old, xmax_old, ymin_old, ymax_old] = axis()
        if xmin == None:
            xmin = xmin_old        
        if xmax == None:
            xmax = xmax_old
        if ymin == None:
            ymin = ymin_old
        if ymax == None:
            ymax = ymax_old
        self.ax.set_xlim( xmin, xmax )
        self.ax.set_ylim( ymin, ymax )
    
    def scale_axis_limits( self, k1, k2, k3, k4 ):
        # k1,   k2,   k3,   k4 are scaling parameters for
        # xmin, xmax, ymin, ymax, respectively.
        sca( self.ax )
        [xmin, xmax, ymin, ymax] = axis()
        self.ax.set_xlim( xmin * k1, xmax * k2 )
        self.ax.set_ylim( ymin * k3, ymax * k4 )
    
    def show_figure( self ):
        ''' note: if you run this, you won't be able to do anything more with it '''
        #TODO: is the above description correct in these scripts? It certainly appears to be the case when working in the shell.
        self.fig.show()
    
    def save_to_file( self, file_path ):
        self.fig.savefig( file_path, facecolor = self.fig.get_facecolor(), edgecolor = 'none' )
    
    def add_note( self, text, x_location, y_location ):
        self.ax.text( x_location, y_location, text, fontsize = 10, ha = 'center' )
    """
    
    def plot_and_save_models( self, show_plot = False, ):
            '''
            TODO: then, make separate methods for 
                ( 1 ) viewing and 
                ( 2 ) saving to file the figures.
            '''
            # Get data for each model and plot it.
            for golf_ball_model, color, angle in zip( self.golf_model_objects, colors, self.initial_angles_degrees ):
                # gather data from the model object
                xcoords = [vector.x for vector in golf_ball_model.x]
                ycoords = [vector.y for vector in golf_ball_model.x]
                dt = golf_ball_model.dt
                # Draw the trajectory
    #            color = colors[ k % len( colors ) ]
    #            '''NB: we use a modulus operation in the previous
    #            line in case we need to re-use colors. Obviously, this 
    #            ideally would be avoided. But, we don't want a silly 
    #            error like that stopping the entire process.'''
                plot( xcoords, ycoords, linewidth = linewidth, color = color )
                # Put a label on the plot which shows dt with the same color
                # as the plot line
                #TODO: change the order that these labels appear.
                figtext( .6, text_position + 0.05, "dt = %s, angle = %s" % ( dt, angle ), fontsize = 10, ha = 'left', color = color )
    #            text( text_size, text_position, "dt = %s, angle = %s" % ( dt, angle ), fontsize = 10, ha = 'left', color = color )
                text_position += text_position_adjustment
                method = golf_ball_model.name()
            # Put a few finishing touches on the plot
            # Display the results
            if show_plot:
                fig.show()
            for file_type in file_types:
"""
