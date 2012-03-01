'''
Created on Feb 7, 2012

@author: res

Assignment:
Write a program to calculate the trajectory of a golf ball and 
calculate the trajectories as a function of angle 
    (use theta = 45, 30, 15 and 9 degrees). 
Choose the initial velocity of the golf ball to be 70 m/s. 
For the drag, assume a general form of
    Fdrag = -C rho A v**2
where 
    rho is the density of air (at sea level), 1.29 kg/m3 , 
    A is the frontal area of the golf ball, 0.0014 m2 , and 
    C is a coefficient to be discussed below.
The following properties have the following impacts:
    dimples -> C is not constant (decreases as v increases).
    spin -> magnus force

For each angle, calculate and compare the trajectories for 
the following cases:
    a) ideal trajectory: no drag and no spin [2 points]
    b) smooth golf ball with drag: choose C = 1/2 [2 points]
    c) dimpled golf ball with drag: 
        choose 
            C = 1/2 for speeds up to v = 14 m/s and 
            C = 7.0/v at higher velocities. 
        This transition to a reduced drag coefficient is due 
        to turbulent flow, caused by the dimples. [3 points]
    d) dimpled golf ball with drag and spin: 
        use a backspin with 
            S0*omega/m = 0.25 s**-1 
        for a typical case. [3 points]
'''

from homework3_Golf import Golf

if __name__ == '__main__':
    file_directory = "/home/res/Documents/duke/2012S/PHY260/homeworks/homework3/"
#    TODO: edit the following line, implement differently when processing data for different parts of the assignment.
    file_name_base = "PHY260_Grisaitis_homework3_"
    file_type_extensions = ['png']
    
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
        
    
