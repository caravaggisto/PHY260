'''
Created on Mar 3, 2012

@author: res

Warning: this code works, I think, but it's pretty crazy to follow.

For better samples, please see the BiotSavart module from this 
assignment, or hw6 or the final. 
'''

# set to True to print diagnostics to console during execution:
DEBUG = False 

import numpy as np
from np import pi
from matplotlib import cm
from matplotlib.axes import rcParams
from matplotlib.pyplot import contour, title, colorbar, show, clabel, plot, \
    colors, figure, xlabel, ylabel, imshow
from Grapher_v1 import Grapher
import datetime
import random

class Dipole_Grid( object ):
    '''
    This class implements the equation describing an electric dipole.
    '''
    def __init__( self, a, R, i, charges ):
        self.a, self.R, self.i = a, R, i
        self.h = float( a ) / i     # space between grid points
        
        n1 = 32     # n for the least dense grid (i=1)
        self.n = ( n1 + 1 ) * i - 1
        # note on `n`:
        # we have to be careful with choosing n, because
        #     1. we need the point charges to fit exactly on the grid
        #     2. we need the boundary condition to hold, namely that 
        #        V(R) = 0 for R > 10.
        n = self.n
        grid = np.zeros( ( n + 2, n + 2 ) )   # make a nice grid matrix, with room for boundaries.
        self.phi = grid.copy()  # charge density, initially zero throughout
        self.rho = grid.copy()  # potential
        # add the point charges to rho:
        self.rho[( n + 1 ) / 2][i * n1 / 2] = charges[0] # integer division done on purpose
        self.rho[( n + 1 ) / 2][i * ( n1 / 2 + 1 )] = charges[1]

        self.inner_points = self.get_inner_pts_bool()
        self.data_about_runs = []
        #TODO: (not urgent) develop a new system of organizing results of models...
        # probably different subclass structures...
    
    def distance_from_origin( self, i, j ):
        h, n = self.h, self.n
        r = h * np.sqrt( ( i - 0.5 * n ) ** 2 + ( j - 0.5 * n ) ** 2 )
        return r
    
    def get_V_of_r( self ):
        ''' returns array with (r,V(r)) coordinates '''
        r, V = [], []
        pts_to_consider = self.get_list_of_pts_inside_boundary()
        for ( i, j ) in pts_to_consider:
            r.append( self.distance_from_origin( i, j ) )
            V.append( self.phi[i][j] )
        return r, V
    
    def is_inside_boundary( self, i, j ):
        ''' Boolean, if grid[i][j] is inside the boundary condition '''
        return self.distance_from_origin( i, j ) < self.R
    
    def get_list_of_pts_inside_boundary( self ):
        # the boy scout way of doing it:
        pts = []
        for i in range( 1, self.n + 1 ):
            for j in range( 1, self.n + 1 ):
                if self.is_inside_boundary( i, j ):
                    pts.append( ( i, j ) )
        return pts
    
    def get_inner_pts_bool( self ):
        # fancy / quick and dirty way of doing it, using the built-in 
        # numpy function.
        return np.fromfunction( self.is_inside_boundary, self.phi.shape )
    
    def jacobi_helper( self, i, j ):
        phi_old, h = self.phi, self.h
        return 0.25 * ( phi_old[i + 1][j] \
                      + phi_old[i - 1][j] \
                      + phi_old[i][j + 1] \
                      + phi_old[i][j - 1] \
                      + h ** 2 * self.rho[i][j] )

    def solve_jacobi( self, epsilon, N_iter_max ):
#        h, n = self.h, self.n
        n = self.n
        N_iter = 0
        self.N_iters = []
        self.errors = []
        error = epsilon + 1     # this is just so we can enter the while loop below:
        #TODO: ramp up N_iter limit for better performance
        while error > epsilon and N_iter < N_iter_max:
            phi_old = self.phi.copy()
            error = 0
            for i in range( 1, n + 1 ):
                for j in range( 1, n + 1 ):
                    if self.inner_points[j][i]:
                        self.phi[i][j] = self.jacobi_helper( i, j )
                        error += abs( self.phi[i][j] - phi_old[i][j] )
            N_iter += 1
            error /= float( n + 1 ) ** 2
            # TODO: comment these out when doing part C
            self.N_iters.append( N_iter )
            self.errors.append( error )
            if DEBUG:
                print "N_iter, error = %d, %.1e" % ( N_iter, error )
        self.N_iter = N_iter
        self.epsilon = epsilon
        return self.N_iters, self.errors
    
    def accuracy( self, phi_old ):
        ''' for SOR, as long as this value is lower than 0.1,
        I'm able to relax the density at which I run the model. '''
        phi_new = self.phi
        fractional_changes = np.divide( phi_old, phi_new ) - 1
        # n.b.: addition (or subtraction) is performed element-wise 
        # by numpy arrays.
        biggest = np.absolute( fractional_changes ).max()
        return biggest
#        return np.abs( ( v_new - v_old ) / v_new )

    def less_dense_grid( self, k ):
        n = self.n
        grid = self.phi.copy()
        grid = grid.astype( bool )
        grid.fill( False )
        offset = random.randint( 0, k - 1 )
        for i in range( 1, n + 1 ):
            for j in range( 1, n + 1 ):
#                if ( ( i % k ) + ( j % k ) ) == 0:
                if ( ( ( i + offset ) % k ) + ( ( j + offset ) % k ) ) == 0:
                # integer division done on purpose...
                # This cuts number of sites from n^2 
                # to (n/k)^2
                    grid[i][j] = True
        return grid
    
    def solve_SOR( self, epsilon, N_iter_max, accuracy_threshold, w = None ):
        ''' Successive Over-Relaxation: a modification to Jacobi Relaxation '''
        n, N_iter = self.n, 0
        if w == None:
            # if no w given, then use the optimal value:
            w = 2.0 / ( 1 + pi / ( n ) ** 0.5 )
        self.N_iters, self.errors = [], []
        error = epsilon + 1     # so we can enter the while loop below:
        while error > epsilon and N_iter < N_iter_max:
            k = 1 # initial density
            grid_SOR = self.less_dense_grid( k )
            phi_old = self.phi.copy()
            error = 0
            for i in range( 1, n + 1 ):
                for j in range( 1, n + 1 ):
                    if self.inner_points[i][j] and grid_SOR[i][j]:
                        self.phi[i][j] = ( 1 - w ) * self.phi[i][j] \
                                        + w * self.jacobi_helper( i, j )
                        error += abs( self.phi[i][j] - phi_old[i][j] )
            N_iter += 1
            error /= float( n + 1 ) ** 2
            if self.accuracy( phi_old ) < accuracy_threshold:
                k += 1
            elif k > 1:
                k -= 1
            self.N_iters.append( N_iter )
            self.errors.append( error )
            if DEBUG:
                print "N_iter, error = %d, %.2e" % ( N_iter, error )
                print error > epsilon
                print N_iter < N_iter_max
        self.N_iter = N_iter
        self.epsilon = epsilon
        return self.N_iters, self.errors
    
    def save_data_to_file( self, file_path = None ):
        # hmm... perhaps it'd be useful to use a local SQLite database here.
        if file_path == None:
            file_path = '/tmp/phy260_midterm-i%d-N_iter%d-ep%d.txt' % ( self.i, self.N_iter, self.epsilon )
        np.save( file_path, self.phi )
    
    def load_data_from_file( self, file_path = '/tmp/tmp.txt.npy' ):
        self.phi = np.load( file_path )
    
    def plot_results( self ):
        h, n = self.h, self.n
        # next few lines auto-generate arrays for plotting
        grid_side = np.linspace( -( n - 1 ) * h / 2, ( n - 1 ) * h / 2, n + 2 )
        X, Y = np.meshgrid( grid_side, grid_side )
        # equipotential values to show in the figure:
        levels = np.array( [ -1e-4, -1e-5, -1e-7, -1e-10, 1e-10, 1e-7, 1e-5, 1e-4] )
        # initialize a figure, housekeeping.
        fig = figure( figsize = ( 7.5, 7.5 ) )
        fig.patch.set_facecolor( 'white' )
        # add a gradient background for potential:
        im = imshow( self.phi, interpolation = 'bilinear',
                cmap = cm.RdBu, extent = ( -10, 10, -10, 10 ) )
        # add equipotential lines
        CS = contour( X, Y, self.phi, levels,
            # `colors` cycles from red to blue to represent equipotential lines
            colors = [cm.RdBu( k * 1.0 / 9 ) for k in range( 1, 9 )],
            linewidths = 4,
            )
        # add title, axis labels:
        title( 'Electrostatic Potential, dipole. Method: SOR.\nN_iter=%d, n=%d, epsilon=%1.e' \
         % ( self.N_iter, self.n, self.epsilon ) )
        xlabel( 'distance' )
        ylabel( 'distance' )
        rcParams['font.size'] = 10
        # make two color bars, to explain the figure
        colorbar( CS, format = '%.e', extend = 'both', orientation = 'vertical', drawedges = True )
        colorbar( im, format = '%.e', extend = 'both' , orientation = 'horizontal' )
        clabel( CS, fmt = '%.e',
#                manual = True, 
                inline_spacing = 0 )
        return fig
        
    def plot_V_of_r( self ):
        r, V = self.get_V_of_r()
        graph = Grapher( "V(r): Electric Potential of a Static Dipole", '', 'radius from origin', 'Potential', 'upper right' )
        graph.add_data_scatter( r, V, 'V(r)' )
        graph.save_to_file( '/home/res/Documents/duke/2012S/PHY260/midterm/V_r.png' )
        graph.show_figure()
        
    def plot_N_iter_with_error( self ):
        N_iters, errors = self.N_iters, self.errors
        N_iter_with_error = Grapher( "V(r): N_iter vs. error", '', 'N_iter, the number of relaxations', 'error', 'upper right' )
        N_iter_with_error.add_data_scatter( N_iters, errors, 'Error = Sum(squared differences of phi[i][j], between iterations' )
#        N_iter_with_error.ax.set_ylim( 0, 2e-5 )
        N_iter_with_error.save_to_file( '/tmp/N_iter_with_error2.png' )
        N_iter_with_error.show_figure()

##################################################################
# Generate results!
##################################################################


##############
# Constants for all the models
##############
a = .6
R = 9.8     # some parts of the grid are barely under 10 points from the origin...
i = 1
charges = ( -1, 1 )
epsilon = 1e-9
N_iter_max = 100000
# for use with SOR:
accuracy = 0.1

##############
# get time for file labeling, whenever we want it
##############
time_string = lambda : datetime.datetime.now().strftime( '%m-%d-%Y_%H-%M-%S' )
# format: '12-25-1995_11-38-05'
"""
##############
### Do Jacobi stuff
##############
grid1 = Dipole_Grid( a, R, i, charges )
jac_N_iters, jac_errors = grid1.solve_jacobi( epsilon, N_iter_max )
#grid1.save_data_to_file( '/tmp/jacobi_' + time_string() )
#grid1.load_data_from_file()
############################## A A A A A A A A A##############
#fig1 = grid1.plot_results()
#fig1.show()

############################## B B B B B B B B B##############
N_iters, errors = jac_N_iters, jac_errors
graph = Grapher( "V(r): N_iter vs. error", '', 'N_iter, the number of relaxations', 'error', 'upper left' )
graph.add_data_scatter( N_iters, errors, 'Sum of squared errors' )
#        N_iter_with_error.ax.set_ylim( 0, 2e-5 )
#graph.add_data( N_iters, errors, 'Sum of squared errors' )
graph.save_to_file( '/tmp/jacobi_' + time_string() + '.png' )
#graph.show_figure()

############################## C C C C C C C C C##############
# Compare results:
#     N_iter with n

##############
### Do SOR stuff
##############
grid2 = Dipole_Grid( a, R, i, charges )
sor_N_iters, sor_errors = grid2.solve_SOR( epsilon, N_iter_max, accuracy )
if DEBUG:
    print "-----------------------------------------"
    print grid2.N_iter
    print "-----------------------------------------"
#grid2.save_data_to_file()
#grid2.load_data_from_file()
fig2 = grid2.plot_results()
fig2.show()
#grid2.plot_V_of_r()
#grid2.plot_N_iter_with_error()
graph = Grapher( "V(r): N_iter vs. error", '', 'N_iter, the number of relaxations', 'error', 'upper left' )
graph.add_data_scatter( sor_N_iters, sor_errors, 'Sum of squared errors' )
#        N_iter_with_error.ax.set_ylim( 0, 2e-5 )
#graph.add_data( N_iters, errors, 'Sum of squared errors' )
graph.save_to_file( '/tmp/jacobi_' + time_string() + '.png' )
graph.show_figure()
"""
##############
# make set of i's for which to evaluate
##############
n_from_i = lambda i: 33 * i - 1
i_list = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 )
n_list = map( n_from_i, i_list )

# make grids
jac_grids = []
sor_grids = []
for i in i_list:
    jac_grids.append( Dipole_Grid( a, R, i, charges ) )
    sor_grids.append( Dipole_Grid( a, R, i, charges ) )

jac_N_iters_list = []
sor_N_iters_list = []
# run models
for jac_grid in jac_grids:
    jac_grid.solve_jacobi( epsilon, N_iter_max )
#    print 'hi'
    jac_N_iters_list.append( jac_grid.N_iter )

    
for sor_grid in sor_grids:
    sor_grid.solve_SOR( epsilon, N_iter_max, accuracy )
#    print 'hey'
    sor_N_iters_list.append( sor_grid.N_iter )

# make plot figure of results
title = 'N iterations vs. n | SOR, Jacobi compared.'
subtitle = ''
x_label = 'n = sites on grid'
y_label = 'N = number of iterations'
legend_loc = 'upper left'
graph = Grapher( title, subtitle, x_label, y_label, legend_loc )
graph.add_data( n_list, jac_N_iters_list, 'Jacobi' )
graph.add_data( n_list, sor_N_iters_list, 'SOR' )
graph.save_to_file( '/home/res/Documents/duke/2012S/PHY260/midterm/jac_sor_' + time_string() + '.png' )

# save the data to be plotted
np.save( '/home/res/Documents/duke/2012S/PHY260/midterm/jac_N_iter_vs_n_' + time_string(), np.asarray( zip( n_list, jac_N_iters_list ) ) )
np.save( '/home/res/Documents/duke/2012S/PHY260/midterm/sor_N_iter_vs_n_' + time_string(), np.asarray( zip( n_list, sor_N_iters_list ) ) )

# show the plot
graph.show_figure()

print "Done with midterm_Dipole.py!"
