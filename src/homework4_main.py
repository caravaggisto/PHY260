'''
Created on Feb 7, 2012

@author: res



Assignment: Oscillatory Motion and Chaos

This is a test I wonder if this wrapper things really works and if obama will get relected and if my artichoekes will taste as good as I hope and if Elena really likes me ;)

'''

from homework4_Pendulum import Pendulum    
import numpy
from numpy import log10, pi, linspace, cos, sin
import ODE_Solver_v3
from matplotlib.pyplot import figure, plot, gcf, show

#file_directory = "/home/res/Documents/duke/2012S/PHY260/homeworks/homework4/"
#file_name_base = "PHY260_Grisaitis_homework3_"
dir_to_save = "/tmp/"
name_base = "test_"
file_type = 'png'

'''
# Test case: u = cos(t)
import ODE_Solver_v3
from matplotlib.pyplot import figure, plot, savefig
U0, omega_D, alpha_D, l, gamma, linear = [0, 0], 0.666, 0.2, 9.8, 0.25, True
#TODO: add more driving frequencies later. First, just make it work.
f = Pendulum( U0, omega_D, alpha_D, l, gamma, linear )
u_init = U0    # initial condition
nperiods = 10     # no of oscillation periods
T = 2 * pi * nperiods    # update this to work with omega_D, and eventually many values of omega_D to be graphed together.
for method_class in ODE_Solver_v3.ForwardEuler, ODE_Solver_v3.RungeKutta2:
    if method_class == ODE_Solver_v3.ForwardEuler:
        npoints_per_period = 200
    elif method_class == ODE_Solver_v3.RungeKutta2:
        npoints_per_period = 20
    n = npoints_per_period * nperiods
    t_points = linspace( 0, T, n + 1 )      # should I make this a numpy array?
    method = method_class( f )
    method.set_initial_condition( U0 )
    u, t = method.solve( t_points )
    # u(t) is a 2 x n array with [u0,u1] for all t's
    
    # get the u0 values from u for plotting
    u0_values = u[:, 0]
    u1_values = u[:, 1]
    u0_exact = cos( t )
    u1_exact = -sin( t )
    fig_u0 = figure()
    alg = method_class.__name__  # (class) name of algorithm
    plot( t, u0_values, 'r-',
         t, u0_exact, 'b-',
#         legend = ( 'numerical', 'exact' ),
#         title = 'Oscillating system; position - %s' % alg )
 )
    fig_u0.savefig( dir_to_save + name_base + '%s.' % alg + file_type )
    fig_u1 = figure()
    plot( t, u1_values, 'r-',
         t, u1_exact, 'b-',
#         legend = ( 'numerical', 'exact' ),
#         title = 'Pendulum; angle - %s' % alg )
 )
    fig_u1.savefig( dir_to_save + name_base + '%s.' % alg + file_type )
'''




''' Part A 
in document
'''
''' Part B 
For 10 driving frequencies, including resonance:
    Plot 1: functions of time
    - theta(t), omega(t)
        * via Euler-Cromer
        * via Runge-Kutta 2nd order
    Plot 2: functions of the driving freq
    - amplitude
    - phase shift
    - full-width at half-maximum
        * compare this with gamma (damping coefficient)
'''

#TODO: after debugging, decrease the time step for better resolution

#defaults:
( U0, __dummy__, alpha_D, length, gamma, linear ) = \
( [0.2, 0], 0.667, 0.5, 9.8, 0.25, True )

#TODO: uncomment the following line when stuff works...
#driving_frequencies = numpy.logspace( log10( 0.666 ) - 3, log10( 0.666 ) + 6, 10 )
# = array([0.000666,..., 0.666,..., 666000])
#driving_frequencies = numpy.logspace( -1, 1, 3 )
driving_frequencies = [0.667]
# = array([ 0.1, 1, 10 ])

pendula = []
for omega_D in driving_frequencies:
    for U0 in [[0.2, 0], [-3.14, 0]]:
        pendula.append( Pendulum( U0, omega_D, alpha_D, length, gamma, linear ) )

nperiods = 10     # no of oscillation periods
T = 2 * pi * nperiods    # update this to work with omega_D, and eventually many values of omega_D to be graphed together.
#TODO: remove the following after debugging
T = 60

for method_class in ODE_Solver_v3.ForwardEuler, ODE_Solver_v3.RungeKutta2:
    npoints_per_period = 200
    n = npoints_per_period * nperiods
    t_points = linspace( 0, T, n + 1 )      # should I make this a numpy array?
#    figure( 'theta' )       # make figure for theta(t)
    fig_u0 = figure( 1 )       # make figure for theta(t)
#    figure( 'omega' )       # and for omega(t)
    fig_u1 = figure( 2 )       # make figure for theta(t)
    #plot models for each driving frequency
    for f in pendula:
        method = method_class( f )
        method.set_initial_condition( f.U0 )
        #TODO: how do I want to store results for comparing results by omega_D?
        u, t = method.solve( t_points )
        # u(t) is a 2 x n array with [u0,u1] for all t's
        u0_values = u[:, 0]  # get the u0 values from u for plotting
        u1_values = u[:, 1]
        figure( fig_u0.number )
        plot( t, u0_values, 'r-' )
        figure( fig_u1.number )
        plot( t, u1_values, 'r-' )
    figure( fig_u0.number )
    u0_exact = cos( t )
    # plot cos(t) as blue:
#    plot( t, u0_exact, 'b-' )
    figure( fig_u1.number )
    u1_exact = -sin( t )
#    plot( t, u1_exact, 'b-' )
    alg = method_class.__name__  # (class) name of algorithm
#    plot( t, u0_values, 'r-',
#         t, u0_exact, 'b-',
##         legend = ( 'numerical', 'exact' ),
##         title = 'Oscillating system; position - %s' % alg )
# )
    #TODO: edit title, legend, notes for each figure.
#         legend = ( 'numerical', 'exact' ),
#         title = 'Pendulum; angle - %s' % alg )
    for fig in fig_u0, fig_u1:
        fig.savefig( dir_to_save + name_base + str( fig.number ) + '%s.' % alg + file_type )

''' Part C 
For omega = resonance frequency:
    Plot 1: functions of time
    - Energy(t)
        * potential
        * kinetic
        * total
        *** Note: over 10 periods (= 10*2pi/omega)
'''
''' Part D 
For omega = resonance frequency:
    Plot 1: functions of time
    - get from part B: theta(t), omega(t) w/ linear restoring, alpha_D = 0.2 (default)
    - theta(t), omega(t) w/ NON-linear restoring, alpha_D = 0.2 (default)
    - theta(t), omega(t) w/ NON-linear restoring, alpha_D = 1.2 (increased)
        * via Euler-Cromer
        * via Runge-Kutta 2nd order
'''
''' Part E 
For omega = 0.666 s-1:
    For alpha_D = 0.2, 0.5 and 1.2:
        For initial angles slightly off-center (dtheta~.001rad):
            Compute the difference in theta between the models as a function of time
            Plot the results and estimate the Lyapunov exponent of the system.
'''




if __name__ == '__main__':
    file_directory = "/home/res/Documents/duke/2012S/PHY260/homeworks/homework3/"
#    TODO: edit the following line, implement differently when processing data for different parts of the assignment.
    file_name_base = "PHY260_Grisaitis_homework3_"
    file_type_extensions = ['png']
    
    # Parameters to examine:
    initial_angles_degrees = ( 9, 15, 30, 45 )
    initial_speed = 70
    dt = 0.01
    
    # This array contains the parameters specific to each part of the problem.
    question_parameters = [ \
        ['A', 0.0, False, False, 'Ideal'], \
        ['B', 0.5, False, False, 'With air resistance'], \
        ['C', 0.5, True, False, 'With air resistance and dimples'], \
        ['D', 0.5, True, True, 'With air resistance, dimples, and spin']]
    
    for ( part, C, dimpled, spin, info ) in question_parameters:
        # initialize a model object for each part of the problem
        golf_ball = Golf( initial_angles_degrees, initial_speed, C, dimpled, spin )
        golf_ball.solve_with_Euler( dt )
        file_name_base = "PHY260_Grisaitis_homework3_" + part
        # plot the model and save results to file
        golf_ball.plot_and_save_models( \
            file_directory = file_directory, \
            file_name_base = file_name_base, \
            file_types = file_type_extensions, \
            info_text = ( ( 'Part %s:' % part ) + info ) \
            )
        
    
