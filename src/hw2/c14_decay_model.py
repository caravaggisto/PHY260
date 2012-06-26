'''
Created on Jan 31, 2012

@author: res
'''
from PHY260_Grisaitis_homework2_Radioactive_Decay import Radioactive_Decay

class C14():
    '''
    initializes some Carbon-14. We can do all sorts of stuff to it, like:
        simulate_radioactive_decay(initial_moles)
        
    '''
    def __init__( self, initial_mass ):
        # physical contants
        self.c14_molar_mass = 14.003241  # grams per mole (source: http://en.wikipedia.org/wiki/Carbon-14 )
        self.c14_half_life = 5700        # years
        self.c14
    
    def simulate_radioactive_decay( self, initial_moles ):
        decay_simulation = Radioactive_Decay( self.c14_half_life, )
