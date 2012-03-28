'''
Created on Feb 7, 2012

@author: res



Assignment: Oscillatory Motion and Chaos

This is a test I wonder if this wrapper things really works and if obama will get relected and if my artichoekes will taste as good as I hope and if Elena really likes me ;)

'''

from homework4_Pendulum import Pendulum    
import numpy
from numpy import log2, pi, linspace, cos, sin
import ODE_Solver_v3
from matplotlib.pyplot import figure, plot, gcf, show, close, legend, figlegend, \
    gca
from Grapher_v1 import Grapher

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
U0 = [0.2, 0]
driving_frequencies = [0.667]
#driving_frequencies = numpy.logspace( log2( 0.666 ) - 3, log2( 0.666 ) + 6, 10, base = 2 )
# = array([0.000666,..., 0.666,..., 666000])
alpha_D = 1.2
length = 9.8
gamma = 0.25
linear = False

# construct pendulum objects
pendula = []
for omega_D in driving_frequencies:
    pendula.append( Pendulum( omega_D, alpha_D, length, gamma, linear ) )

nperiods = 10     # no of oscillation periods
T = 2 * pi * nperiods    # update this to work with omega_D, and eventually many values of omega_D to be graphed together.
#TODO: remove the following after debugging
T = 100

for method_class in ODE_Solver_v3.EulerCromer, ODE_Solver_v3.RungeKutta2:
    npoints_per_period = 500
    n = npoints_per_period * nperiods
    t_points = linspace( 0, T, n + 1 )      # should I make this a numpy array?
    u0_graph = Grapher( 'Theta(t)', '', 't (seconds)', 'theta (radians)', 'bottom right' )
    #title, subtitle, x_label, y_label, legend_loc 
    u1_graph = Grapher( 'Omega(t)', '', 't (seconds)', 'omega (radians per second)', 'bottom right' )
    for f in pendula:
        method = method_class( f )
        method.set_initial_condition( U0 )
        #TODO: how do I want to store results for comparing results by omega_D?
        u, t = method.solve( t_points )
        # u(t) is a 2 x n array with [u0,u1] for all t's
        u0_values = u[:, 0]  # get the u0 values from u for plotting
        u0_graph.add_data( t, u0_values, str( f.omega_D ) )
        u1_values = u[:, 1]
        u0_graph.add_data( t, u1_values, str( f.omega_D ) )
#    u0_graph.show_figure()
    u1_graph.show_figure()

for method_class, color in [( ODE_Solver_v3.RungeKutta2, 'r-' ), ( ODE_Solver_v3.EulerCromer, 'b-' )]:
    npoints_per_period = 500
    n = npoints_per_period * nperiods
    t_points = linspace( 0, T, n + 1 )      # should I make this a numpy array?
#    figure( 'theta' )       # make figure for theta(t)
    fig_u0 = figure( 1 )       # make figure for theta(t)
#    figure( 'omega' )       # and for omega(t)
    fig_u1 = figure( 2 )       # make figure for theta(t)
    #plot models for each driving frequency
    for f in pendula:
        method = method_class( f )
        method.set_initial_condition( U0 )
        #TODO: how do I want to store results for comparing results by omega_D?
        u, t = method.solve( t_points )
        # u(t) is a 2 x n array with [u0,u1] for all t's
        u0_values = u[:, 0]  # get the u0 values from u for plotting
        u1_values = u[:, 1]
        #u0_graph.add_data(t, u, legend_label)
        #u1_graph.add_data(t, u, legend_label)
        figure( fig_u0.number )
        plot( t, u0_values, color )
        figure( fig_u1.number )
        plot( t, u1_values, color )
    figure( fig_u0.number )
    fig_u0.get_axes()[0]
    legend()
    figure( fig_u1.number )
    alg = method_class.__name__  # (class) name of algorithm
#    plot( t, u0_values, 'r-',
#         t, u0_exact, 'b-',
##         legend = ( 'numerical', 'exact' ),
##         title = 'Oscillating system; position - %s' % alg )
# )
#         legend = ( 'numerical', 'exact' ),
#         title = 'Pendulum; angle - %s' % alg )
    for fig in fig_u0, fig_u1:
        fig
        fig.savefig( dir_to_save + name_base + str( fig.number ) + '%s.' % alg + file_type )
print gcf().__hash__
show()


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

