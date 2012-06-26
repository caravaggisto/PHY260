'''
PHY260 | homework 2
due: Jan 30, 2011

@author: William C Grisaitis
'''

from PHY260_WilliamGrisaitis_homework2_Radioactive_Decay import Radioactive_Decay

if __name__ == '__main__':
    Carbon14_half_life = 5700 # years
#    TODO update molar mass, be exact, for Carbon-14
    Carbon14_molar_mass = 14.0 # kg per mol
    initial_mass = 10e-12 # kg
    
    decay = Radioactive_Decay( Carbon14_half_life, initial_mass, Carbon14_molar_mass )
    
    ''' part b: Numerically calculate the activity of the sample 
    over a duration of 20,000 years. Use numerical time-width 
    steps of 10 and 100 years. Plot results in appropriate units 
    along with the analytical (exact) result - all in one plot.
    '''
    
    ''' part c: Re-calculate with time width of 1,000. Re-plot. 
    Is the accuracy still acceptable? What is the percentage 
    deviation after two years from the analytical result? Is the 
    error as large as you would expect, given that we are omitting 
    a second order term?
    '''
    
    pass
