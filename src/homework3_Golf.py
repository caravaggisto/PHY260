'''
Created on Feb 7, 2012

@author: res
'''

from numpy import pi, sin, cos
from Forces import Force_Gravity, Ball_Drag, Magnus_Effect
from Vector_3D import V3
from ODE_Solver import EulerSolver, AboveGround
from matplotlib.pyplot import plot, text, title, xlabel, ylabel, show, figure, \
    xlim, ylim, figtext

class Golf( object ):
    '''
    A class to implement different models of golf ball trajectory.
    '''
    def __init__( self, initial_angles_degrees, initial_speed = 70.0, C = 0, dimpled = False, spin = False ):
        '''
        Needed:
            initial_angles_degrees
        Optional (= default):
            initial_speed = 70.0, 
            C = 0, 
            dimpled = False, 
            spin = False
        '''
        # Initial conditions
        self.initial_angles_degrees = initial_angles_degrees
        self.initial_angles_in_radians = [( angle * pi / 180.0 ) for angle in initial_angles_degrees ]
        self.x0 = V3( 0.0, \
                      0.0, \
                      0.0 )
        self.initial_speed = initial_speed
        self.v0 = None
        self.C = C              # the drag coefficient.
        self.dimpled = dimpled  # a boolean: Is the ball dimpled?
        self.spin = spin        # a boolean: Does the ball have spin?
        # a few constants
        self.mass = 0.0459      # mass of golf ball, kg
        self.g = 9.81           # gravity, m / s**2
        self.rho = 1.29         # density of air at sea level, kg / m**3
        self.A = 0.0014         # frontal area of the golf ball, m**2
        self.S0_times_omega_over_m = 0.25
        # an empty array to store model results
        self.golf_model_objects = []
        self.results = []
        self.plots = []
    
    def terminate( self, u, t, step_no ):
        return False if u[step_no, 2] >= 0 else True
    
    ''' Define Forces '''
    
    def force_gravity( self ):
        # Force of gravity
        g, m = self.g, self.mass
        return Force_Gravity( g, m )
    
    def force_air_resistance( self ):
        # Force of air resistance
        A, rho, C, dimpled = self.A, self.rho, self.C, self.dimpled
        return Ball_Drag( A, rho, C, dimpled )
    
    def force_magnus( self ):
        # Force due to Magnus Effect
#        S0 * omega / mass = 0.25
        S0_times_omega_over_m = self.S0_times_omega_over_m
        S0_times_omega = S0_times_omega_over_m * self.mass
        if self.spin:
            omega_vector = V3( 0.0, 0.0, 1.0 )
        else:
            omega_vector = V3( 0.0, 0.0, 0.0 )
        return Magnus_Effect( S0_times_omega, omega_vector )

#        return Magnus_Effect( S0_times_omega_over_m, omega_vector )
    
    def solve_with_Euler( self, dt ):
        m = self.mass
        Fg = self.force_gravity()
        Far = self.force_air_resistance()
        Fm = self.force_magnus()
        # Run the simulation
        for angle in self.initial_angles_in_radians:
            golf_model = EulerSolver( Fg + Far + Fm, m )
            self.golf_model_objects.append( golf_model )
            self.v0 = V3( self.initial_speed * cos( angle ), \
                          self.initial_speed * sin( angle ), \
                          0.0 )
            golf_model.run( self.x0, self.v0, dt, AboveGround() )
            # Extract the trajectory in the XY plane
            xcoords = [vector.x for vector in golf_model.x]
            ycoords = [vector.y for vector in golf_model.x]
            # Save the data for this model in the local results array
            self.results.append( [ zip( xcoords, ycoords ), dt, angle ] )
            print '%6.2f' % xcoords[-1]
    
    def plot_and_save_models( self, \
                              linewidth = 1.0, \
                              colors = ['r', 'y', 'g', 'b'], \
                              text_size = 10, \
                              text_position = 0.43, \
                              text_position_adjustment = 0.09, \
                              file_directory = "/tmp/", \
                              file_name_base = "PHY260_test", \
                              file_types = ['png'], \
                              show_plot = False, \
                              info_text = '' ):    
        # Turn on the grid on the plot
        fig = figure( figsize = ( 7.5, 7.5 ) )
        ax = fig.add_subplot( 111 )
        ax.grid( True )
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
        title( 'Simulated Trajectories of a Golf Ball (%s Method)' % method )
        figtext( .5, .85, info_text, fontsize = 10, ha = 'center' )
        xlabel( 'x (m)' )
        ylabel( 'y (m)' )
        xlim( 0, 500 )
        ylim( 0, 150 )
        fig.patch.set_facecolor( 'white' )
        # Display the results
        if show_plot:
            fig.show()
        for file_type in file_types:
            fig.savefig( file_directory + file_name_base + "." + file_type, facecolor = fig.get_facecolor(), edgecolor = 'none' )

    
    
    
