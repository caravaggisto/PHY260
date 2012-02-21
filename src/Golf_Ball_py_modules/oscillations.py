#!/usr/bin/env python

from cpode import *
from cpforces import *
from pylab import *
import math

# Prepare the force acting on the test mass (m = 1.0)
springCoefficient = 1.0
force = RestoringForce(springCoefficient, 0.0)

# Initial conditions
x0 = 0.0
v0 = 1.01

# When do we stop?
oscillationPeriod = 2*math.pi/sqrt(springCoefficient)
timeLimit = 10*oscillationPeriod

# The tuple of ODE algorithms that we are going to use
odeSolvers = (RK4, EulerCromer, EulerSolver)

# Run for a set of different time deltas.
# We will display them using different colors.
# Run "help(colors)" at the prompt to see
# the inctructions on the use of colors.
timedeltas = (0.02,  0.01)
plotcolors = ('r',   'g')

# Cycle over the ODE algorithms
for algorithm in odeSolvers:
    solver = algorithm(force)

    # Make a new plot and turn on the grid on the plot
    ax = figure().add_subplot(111)
    ax.grid(True)

    text_position = 0.03

    # Cycle over time deltas and corresponding colors
    for dt, color in zip(timedeltas, plotcolors):
        # Run the simulation
        solver.run(x0, v0, dt, TimeLimit(timeLimit))

        # Draw the displacement
        plot(solver.t, solver.x, linewidth=1.0, color=color)

        # Put a label on the plot which shows dt with the same color
        # as the plot line
        text(0.05, text_position, "dt = %s" % dt, color=color,
             transform = ax.transAxes)
        text_position += 0.05

    # Put a few finishing touches on the plot for this solver
    title('Simulated Oscillations (%s Method)' % solver.name())
    xlabel('t (s)')
    ylabel('x (m)')

# Display the results
show()
