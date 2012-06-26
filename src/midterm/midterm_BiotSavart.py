'''
Created on Mar 4, 2012

@author: res
'''

import numpy as np
from numpy import pi
from Grapher_v1 import Grapher
import datetime

##############
# get time for file labeling, whenever we want it
##############
time_string = lambda : datetime.datetime.now().strftime( '%m-%d-%Y_%H-%M-%S' )
# format: '12-25-1995_11-38-05'

class Square_Current_Loop( object ):
    '''
    This class simulates a square wire in the xy-plane 
    centered at the origin, creating a magnetic field 
    symmetric about the xz and yz planes,
        B(x,y,z) = B(-x,y,z) = B(x,-y,z) = B(-x,-y,z)
    and antisymmetric about the xy plane:
        B(x,y,z) = -B(x,y,-z)
    where B(x,y,z) is the magnetic field vector at the 
    (x,y,z) position.

    TODO: make iterable! DONE.
    '''
    def __init__( self, side_length, points_per_side ):
        self.l = float( side_length )
        dl = float( side_length ) / points_per_side   # distance between each successive dr' on the current loop.
        # make an array of 4 triplets.
        verticies = np.zeros( ( 4, 3 ) )
        verticies[0] = np.array( ( 0.5, 0.5, 0 ) )  # a square in the xy-plane :)
        verticies[1] = np.array( ( -0.5, 0.5, 0 ) )
        verticies[2] = np.array( ( -0.5, -0.5, 0 ) )
        verticies[3] = np.array( ( 0.5, -0.5, 0 ) )
        # set each point of the current loop:
        self.points = np.zeros( ( 4 * points_per_side, 3 ) )
        # "4" for 4 side per square, "3" for 3 dimensions
        for i in range( 4 ):
            direction = verticies[i] - verticies[i - 1] 
            # = the direction in which we're constructing more points along the side of the current loop square.
            # Note how this is automatically normalized due to how `verticies` is constructed.
            point = verticies[i - 1] * side_length  
            # give the vertex magnitude
            for k in range( points_per_side ):
                point += direction * dl
                # inch along the side of the current loop...
                self.points[ i * points_per_side + k ] = point.copy()
                # go through self.points, changeing 
        self.k = 0 # for iteration
    
    def __iter__( self ):
        return self
    
    def next( self ):
        ''' makes the class iterable for polymorphic goodness' sake '''
        if self.k == len( self.points ):
            raise StopIteration
        else:
            self.k += 1
            return self.points[self.k - 1]
        
class B_Field( object ):
    def __init__( self, mu0_I, loop ):
        self.mu0_I = mu0_I  # a constant
        self.loop = loop    # a Square_Current_Loop object
        
    def db_helper( self, r, this_point, previous_point ):
        ''' db doesn't mean database... rather, a discrete 
        difference in the magnetic field at a point in 
        cartesian space, from 

        this_point      = the coordinate currently being evaluated
        previous_point  = the coord evaluated just previously evaluated.


        '''
        dl_vector = this_point - previous_point
        r_prime = r - this_point
        db = np.cross( dl_vector, r_prime ) * np.vdot( r_prime, r_prime ) ** -1.5
        return db
        
    def __call__( self, r ):
        # where r is a point somewhere in space, at which 
        # we want to know the B(r) vector.
        b = np.zeros_like( r ) # initialize a magnetic field vector
        first_point = self.loop.next() # this is called only once. Thus, 
        # self.loop.next() is, in this case, the first iterable element in the loop object. 
        previous_point = first_point.copy()
        # we need to get the first point before we start iterating,
        # because to get the `dl` vector, we need the previous point
        # visited. Because at the beginning there's no 'previous' point, 
        # we just use the first point twice. This is redundant, and 
        # I should change it, but it's a minor performance issue.
        for this_point in loop:
            b += self.db_helper( r, this_point, previous_point )
            previous_point = this_point.copy()
            # I don't know if this has to be a copy or not, but I'd 
            # rather not take any chances with these schizophrenic numpy objects.
        self.loop.k = 0
        b *= self.mu0_I / ( 4 * pi )
        # slap on the coefficient. 
        return b
        # Et alors, nous avons le Biot-Savart! Hanh hanh hanh...
            
######## 
# Physical constants, given in the assignment:
########
mu0_I = 8 
pi = np.pi # 3.14159265358979323... 
side_length = 1.0

####################################################################
# constants for determining resolution of the model
####################################################################
L = 10.0 # total length in meters of the domain we're interested in.
domain_step = 0.005
axis_points = np.arange( -L * 0.5, L * 0.5 + domain_step, domain_step )
# examine from -L/2 to L/2 in steps of 0.005.

######## 
# Make current loop
########
# make an iterable current loop object that returns 
# successive positions along the current loop.
loop_step = 0.5
points_per_side = int( side_length / loop_step )
loop = Square_Current_Loop( side_length, points_per_side )
b_field = B_Field( mu0_I, loop )

######## 
# Perhaps useless...
list_of_arrays_of_points = []
########

######## 
# Part A: initialize stuff
########
constraints = ( 0, 0, None )
random_sample = axis_points.copy()
random_sample.fill( 0 )
y_values = axis_points.copy()
y_values.fill( 0 )
z_values = axis_points
N = len( z_values ) # = number of points approximated
part_A_points = np.column_stack( ( random_sample, y_values, z_values ) )
list_of_arrays_of_points.append( part_A_points )

######## 
# Part B: initialize stuff
########
constraints = ( None, 0, 1 )
random_sample = axis_points
y_values = axis_points.copy()
y_values.fill( 0 )
z_values = axis_points.copy()
z_values.fill( 1 )
part_B_points = np.column_stack( ( random_sample, y_values, z_values ) )
list_of_arrays_of_points.append( part_B_points )

######## 
# Part C: initialize
########
constraints = ( 0.5, 0, None )
random_sample = axis_points.copy()
random_sample.fill( 0.5 )
y_values = axis_points.copy()
y_values.fill( 0 )
z_values = axis_points
part_C_points = np.column_stack( ( random_sample, y_values, z_values ) )
list_of_arrays_of_points.append( part_C_points )

''' TODO: Optimization. I could combine this entire last chunk of 
 repeated code into a fast-track initializer for any (x,y,z) constraint
 tuple. '''

subtitle = ''
y_label = 'B, magnetic field'
legend_loc = 'upper right'

######## 
# Execute parts of the question
########

# some helper stuff
R = side_length * 0.5 # radius described as half the langth of one side of the current loop.
b_z_circular = lambda z: ( mu0_I * R ** 2 / ( 2 * ( z ** 2 + R ** 2 ) ** ( 1.5 ) ) )

# Part A:
y1_vectors = np.asarray( map( b_field, part_A_points ) )
y2 = map( b_z_circular, axis_points )
for ( label, i ) in zip( ( 'x', 'y', 'z' ), range( 3 ) ):
    title = 'Part A: B(' + label + ') for x = y = 0 (Numerical and Analytical) '
    y1 = np.hsplit( y1_vectors, 3 )[i]   # get z-axis points, in third column
    x_label = label
    graph = Grapher( title, subtitle, x_label, y_label, legend_loc )
    graph.add_data( axis_points, y1, 'Numerical: B for x = y = 0' )
    graph.add_data( axis_points, y2, 'Analytical: B for x = y = 0' )
    graph.save_to_file( '/tmp/' + 'A' + label + '.png' )

# Part B:
y1_vectors = np.asarray( map( b_field, part_B_points ) )
for ( label, i ) in zip( ( 'x', 'y', 'z' ), range( 3 ) ):
    title = 'Part B: B(' + label + ') for y=0, z=1'     
    y1 = np.hsplit( y1_vectors, 3 )[i]   # get x-axis points
    x_label = label
    graph = Grapher( title, subtitle, x_label, y_label, legend_loc )
    graph.add_data( axis_points, y1, 'Approximation: B for y=0, z=1' )
    #graph.show_figure()
    graph.save_to_file( '/tmp/' + 'B' + label + '.png' )

# Part C:
y1_vectors = np.asarray( map( b_field, part_C_points ) )
for ( label, i ) in zip( ( 'x', 'y', 'z' ), range( 3 ) ):
    title = 'Part C: B(' + label + ') for y=0, x=0.5'
    y1 = np.hsplit( y1_vectors, 3 )[i]   # get z-axis points
    x_label = label
    graph = Grapher( title, subtitle, x_label, y_label, legend_loc )
    graph.add_data( axis_points, y1, 'Approximation: B for y=0, x=0.5' )
    graph.save_to_file( '/tmp/' + 'C' + label + '.png' )
    