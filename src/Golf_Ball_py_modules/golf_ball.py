#!/usr/bin/env python

from v3 import *
from cpode import *
from cpforces import *
from pylab import *
import math

# The coordinates follow Giordano & Nakanishi convention:
# x axis is horizontal, along the projectile motion.
# y axis is vertical (directed upwards), and z axis makes
# a left-handed coordinate system with x and y.

# The mass of the golf ball
m = 45.93e-3

# Gravity force
Fg = ForceOfGravity(m)

# Air drag
Fd = GolfBallDrag()

# Magnus force
S0 = 0.25*m
omega = V3(0.0, 0.0, 1.0)
Fm = MagnusForce(S0, omega)

# Prepare the ODE solver
golfBall = EulerSolver(Fg + Fd + Fm, m)

# Initial conditions
x0 = V3(0.0, 0.0, 0.0)
hitAngle = 10*math.pi/180.0
speed = 70.0
v0 = V3(speed*math.cos(hitAngle), speed*math.sin(hitAngle), 0.0)

# Run for a set of different time deltas.
# We will display them using different colors.
# Run "help(colors)" at the prompt to see
# the inctructions on the use of colors.
timedeltas = (1,   0.1,  0.01)
plotcolors = ('g', 'r',   'b')

# Turn on the grid on the plot
ax = figure().add_subplot(111)
ax.grid(True)

text_position = 0.33

# Cycle over time deltas and corresponding colors
for dt, color in zip(timedeltas, plotcolors):
    # Run the simulation
    golfBall.run(x0, v0, dt, AboveGround())

    # Extract the trajectory in the XY plane
    xcoords = [vector.x for vector in golfBall.x]
    ycoords = [vector.y for vector in golfBall.x]

    # Draw the trajectory
    plot(xcoords, ycoords, linewidth=1.0, color=color)

    # Put a label on the plot which shows dt with the same color
    # as the plot line
    text(0.45, text_position, "dt = %s" % dt, color=color,
         transform = ax.transAxes)
    text_position -= 0.09

# Put a few finishing touches on the plot
title('Simulated Trajectories of a Golf Ball (%s Method)' % golfBall.name())
xlabel('x (m)')
ylabel('y (m)')

# Display the results
show()
