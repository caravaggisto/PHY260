'''
Created on Feb 7, 2012

@author: res
'''
from homework3_Golf import Golf

class Problem( object ):
    '''
    classdocs
    '''
    
    def __init__( self, initial_velocity, initial_angles_degrees, filetypes = ['pdf', 'eps', 'png'] ):
        '''
        Constructor
        '''
        self.initial_velocity = initial_velocity
        self.initial_angles_degrees = initial_angles_degrees
        self.filetypes = filetypes
    
    def force_magnus( self, S0wm ):
        pass
    
    def part_a( self ):
        '''
        a) ideal trajectory: no drag (C=0) and no spin (-> Magnus force = 0) [2 points]
        '''
        initial_velocity = self.initial_velocity
        dimpled = False
        spin = 0
        def C( v, dimpled ):
            return 0
        
        
        for initial_angle in self.initial_angles_degrees:
            part_a_ball = Golf( initial_velocity, initial_angle, C, dimpled, spin )
            part_a_ball.display_plot()
            for filetype in self.filetypes:
                part_a_ball.save_plot( filetype )        
    
    def part_b( self ):
        '''
        b) smooth golf ball with drag: choose C = 1/2 [2 points]
        '''
        dimpled = False
    
    def part_c( self ):
        '''
        c) dimpled golf ball with drag: 
            choose 
                C = 1/2 for speeds up to v = 14 m/s and 
                C = 7.0/v at higher velocities. 
            This transition to a reduced drag coefficient is due 
            to turbulent flow, caused by the dimples. [3 points]        '''
        pass
    
    def part_d( self ):
        '''
        d) dimpled golf ball with drag and spin: 
            use a backspin with 
                S0 ω/m = 0.25 s**-1 
            for a typical case. [3 points]        '''
        def dimpled_C ( self, v ):
            '''
            assumes C = 0.5
            '''
            C = self.C
            vmax = 14.0
            if v < vmax:
                return C
            else:
                return vmax / 7
        pass
