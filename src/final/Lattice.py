'''
@author: willg

PHY260, Final Project

This module implements a class structure for a lattice 
for the Ising model.

a) Choose n = 50 and calculate M = N<s> as a function of temperature.
--> create array of temperature values
--> relax metropolis algorithm along that temperature array.
Determine the critical temperature T_C
--> TODO: how to implement finding T_C?

b) Calculate the specific heat per spin (C/N) for the following lattice sizes:
    n in [5,10,20,30,40,50,75,100,200,500]
using the fluctuation-dissipation theorem:
    C = (1/k_B)*(dE/T)^2
and verify that
    Cmax / N ~ log(n)
--> make array of log(n)
--> determine Cmax for each n

    title: "PHY260, Final Exam, A: Total Magnetization as a function of k_B*T.
    x-axis: "k_B*T"
    y-axis: "Total Magnetization = M = N<s>"
    rasterized: False
'''

import numpy as np

class Lattice( object ):
    def __init__( self, J, n, n_steps, kT_array ):
        # Simulation-determined params
        self.J = float( J ) # interaction strength
        self.n = n # for a n x n lattice
        self.n_steps = n_steps
        self.kT_array = kT_array # array of temperatures
        # Constants
        self.dE = np.array( [2, 1, 0, -1, -2] ) * 4
        '''
        # Initialize lattice with random spins
        spins = {0:-1, 1:1}
        lattice_indeces = np.random.randint( 0, 2, ( n, n ) )
        self.lattice = np.zeros( ( n, n ), dtype = int )
        for key, spin in spins.iteritems():
            self.lattice[lattice_indeces == key] = spin
        '''
        # Initialize the lattice with all spin up
        self.lattice = np.ones( ( n, n ), dtype = int )
        self.eTot = 0.0
        
        # Empty arrays to store values for each value of kT
        self.M = np.zeros_like( kT_array ) # to store Magnetizations as fn of kT
        self.C = np.zeros_like( kT_array ) # to store C / N, specific heat capacity
    
    def run_metropolous_algorithm( self ):
        kT_array, n_steps = self.kT_array, self.n_steps
        for temp_index, kT in enumerate( kT_array ):
            print "*"*50, kT
            self.kT = kT
            self.dEB = np.exp( -self.dE / kT )
            print self.dEB
            self.E = []
            self.E_steps = np.zeros( n_steps ) # empty array to store energies for each step
            self.M_steps = np.zeros( n_steps ) # empty array to store magnetizations
            for n_relax in range( n_steps ):
                self.relax( n_relax ) # relaxes the grid `n_steps` times
            self.C[temp_index] = self.calculate_specific_heat_per_spin()
            self.M[temp_index] = self.calculate_M()
    
    def flip( self, n ):
        self.lattice[n] *= -1
    
    def get_neighboring_points( self, coordinates ):
        # returns neighbors of lattice[i,j], assuming 
        # torus-like period boundary conditions
        n = self.n
        i, j = coordinates
        up = ( ( i - 1 ) % n, j )
        down = ( ( i + 1 ) % n, j )
        right = ( i, ( j + 1 ) % n )
        left = ( i, ( j - 1 ) % n )
        return ( up, down, right, left )
    
    def get_dE_index( self, i, j ):
        count = 0
        for n in self.get_neighboring_points( ( i, j ) ):
            # checks if the neighboring spin is the same as [i,j].
            # if so, we increment `count`
            # the higher `count` gets, the lower the energy will be.
            if self.lattice[n] == self.lattice[i, j]:
                count += 1
                # `count` corresponds to the indeces of 
                # self.dE = np.array( [2, 1, 0, -1, -2] )
        return count
    
    def relax( self, n_relax ):
        # flips a random site on the grid.
        n = self.n
        k, l = np.random.randint( 0, n, 2 ) # generates 2 random numbers
        self.flip( ( k, l ) )
        dE_index = self.get_dE_index( k, l )
        # Obtain the according changes in the energy *and* the boltzman stat
        dE, dEB = self.dE[dE_index], self.dEB[dE_index]
        if dE <= 0.0:
            # if there's a decrease in energy, the flip will ALWAYS occur
            # see Giordano, p.214
            if n_relax % 10000 == 0:
                print "%05i: Accepted. dE = %0.2f" % ( n_relax, dE )
            # ... then all we do is add it to the total energy, and continue.
            self.eTot += dE
        else:
            # if energy INcreases... 
            # 1) make a ~Unif(0,1) variable
            x = np.random.random()
            # and if dEB > x, the system successfully increases in energy
            if dEB >= x:
                # 
                # we keep the flip as-is... 
                if n_relax % 10000 == 0:
                    print "%05i: " % n_relax + \
                            "kT: %2.1f, " % self.kT + \
                            "Accepted  Bf: %3.2f, " % dEB + \
                            "%1.4f" % x
                # and add it to the total
                self.eTot += dE
            else:
                # we flip it back
                if n_relax % 10000 == 0:
                    print "%05i: " % n_relax + \
                            "kT: %2.1f, " % self.kT + \
                            "Not Acce. Bf: %3.2f, " % dEB + \
                            "%1.4f" % x
                self.flip( ( k, l ) )
        # Store this energy and magnetization to an array for all relaxations
        self.E_steps[n_relax] = self.eTot
        self.E.append( self.eTot )

    def calculate_M( self ):
        return self.lattice.sum()
    
    def calculate_specific_heat_per_spin( self ):
        ''' C/N ~ [ (J dE) / (kT) ]**2 / n**2
            = Var(E) * ( J / (kT n) )**2 '''
        return np.var( np.asarray( self.E ) ) * ( self.J / ( self.kT * self.n ) ) ** 2 
    
    def get_max_specific_heat( self ):
        return max( self.C )

if __name__ == '__main__':
    '''this code is for debugging the class, making sure 
    the methods are all working the way they're supposed to, etc. '''
    np.random.seed( 1 ) # use the same set of random numbers... forever!
    #===============================================================================
    # Adjustable parameters for the problem 
    J = 1.5
    n = 50 #TODO: make this array-able... in B) we need to compare results for different n.
    n_steps = 20000
    kT_min = 0.0
    kT_max = 6.0
    kT_step = 1.0
    #===============================================================================
    kT_array = np.arange( kT_min, kT_max, kT_step )
    l1 = Lattice( J, n, n_steps, kT_array )
    l1.run_metropolous_algorithm()
    # make lattice data 
    data = {'M vs. T': l1.M}
    # plot the results
    import matplotlib.pyplot as plt
    fig = plt.figure()
    for n, y in enumerate( data ):
        ax = fig.add_subplot( len( data ), 1, n + 1 )
        ax.scatter( l1.kT_array,
                    data[y],
                    c = 'r' )
        ax.set_xlim( kT_min, kT_max )
    plt.savefig( 'test.png' )
