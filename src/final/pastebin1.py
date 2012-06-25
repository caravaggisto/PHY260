import numpy as np

class Lattice_old( object ):
    def __init__( self, J = 1.0, n = 10, kT = 1.0 ):
        self.J = J # interaction strength
        self.n = n # n x n lattice
        self.lattice = np.ones( ( self.n, self.n ), dtype = float )
        self.kT = float( kT )
        self.eTot = 0.0
        
        # The possible values for the change in energy when
        # flipping a spin.
        # All up, the middle one is shifted to down: dE = +8.
        # One do, the middle one is shifted to down: dE = +4.
        # [...]
        # All do, the middle one is shifted to down: dE = -8.
        self.dE = [+8, +4, 0, -4, -8]
        #self.dE = [-8, -4, 0, +4, +8]
        
        self.dEB = np.array( self.dE )
        # Precalculate the corresponding Boltzmann factor
        # for the given temperature.
        self.dEB = np.exp( -self.dEB / self.kT )
        print self.dEB, '****~~~~~******'
        
        # Convenience lists for the properties to be sampled.
        self.E = []
        self.M = []
    
    # Convenience methods.
    def flip( self, i, j ):
        self.lattice[i, j] *= -1
    
    def get_neighboring_points( self, r, c ):
        # r: row, c: column
        # (r, c) is the position of the spin of which
        # the positions of the neighbors are returned.
        # Periodic boundaries:
        # if m == 0, then the neighbor is the element
        # in the last row, i.e. the lattice is
        # on a torus.
        up = [( r - 1 ) % self.n, c]
        le = [r, ( c - 1 ) % self.n]
        ri = [r, ( c + 1 ) % self.n]
        do = [( r + 1 ) % self.n, c]
        return [up, le, ri, do]
    
    # Called upon initialization of the simulation.
    def get_initial_E_total( self ):
        eIni = 0.0
        # Run over the lattice and over neighbors of spin at
        # site 'i,j'.
        for i in range( self.n ):
            for j in range( self.n ):
                for n in self.get_neighboring_points( i, j ):
                    eIni += self.J * self.lattice[i, j] * self.lattice[n[0], n[1]]
        # Devide by 2 to eliminate double counting.
        self.eTot = -eIni / 2.0
    
    def get_dE_index( self, i, j ):
        # Calculating dE depending on number of neighbors which are
        # already down.
        # By convention, initially all spins up and J>0.
        # 'counter' counts how many spins are of like orientation.
        counter = 0
        for n in self.get_neighboring_points( i, j ):
            # We flipped the i,j spin. Now we check its environment.
            # E.g. for all neighbor spins up, and the i,j spin down,
            # 'if' never evaluates to 'True'. We end up with '+8' in
            # dE.
            if self.lattice[n[0], n[1]] == self.lattice[i, j]:
                counter += 1
        # Returns the 'counter' too so the call to the 'index' function
        # can be avoided.
        return counter, self.dE[counter]
    
    def relax( self, n_relax ):
        # Define random site to flip.
        k, l = np.random.randint( 0, self.n - 1 ), np.random.randint( 0, self.n - 1 )
        self.flip( k, l )
        # Obtain the index of the energy value in the prestored
        # energy np.array, and the change in energy for the
        # prospective flip.
        ind, dE = self.get_dE_index( k, l )
        if dE <= 0.0:
            if n_relax % 10000 == 0:
                print "%05i: Accepted. dE = %0.2f" % ( n_relax, dE )
            # Only add change in energy to the total energy.
            self.eTot += dE
        else:
            # Uniform distributed random number.
            x = np.random.random() # x ~ Unif( [0,1) ) 
            # Obtain the corresponding, precalculated Boltzmann
            # factor, avoiding the call to the 'index' builtin
            # function.
            if self.dEB[ind] >= x:
                if n_relax % 10000 == 0:
                    print "%05i: " % n_relax + \
                            "kT: %2.1f, " % self.kT + \
                            "Accepted  Bf: %3.2f, " % self.dEB[ind] + \
                            "%1.4f" % x
                # Adding the change in energy.
                self.eTot += dE
            else:
                # The n_relax is being undone.
                if n_relax % 10000 == 0:
                    print "%05i: " % n_relax + \
                            "kT: %2.1f, " % self.kT + \
                            "Not Acce. Bf: %3.2f, " % self.dEB[ind] + \
                            "%1.4f" % x
                # Set the spin state back.
                self.flip( k, l )
        # Append the total energy of the temperature.
        self.E.append( self.eTot )
        # Append the total magnetization to 'self.M'.  
        self.M.append( np.sum( self.lattice ) )
    
    def heat( self, steps ):
        E = np.array( self.E )
        E2 = E * E
        return 1.0 / self.n * 1.0 / self.n * \
                ( self.J / self.kT ) * ( self.J / self.kT ) * \
                ( 1.0 / steps * np.sum( E2 ) - ( 1.0 / steps * np.sum( E ) ) * ( 1.0 / steps * np.sum( E ) ) )    
    
    def magnet( self, steps ):
        M = np.array( self.M )
        return M / ( steps * self.n ** 2 )
    
    def susz( self, steps ):
        M = np.array( self.M )
        M /= ( self.n * self.n )
        #M /= M[0]
        M2 = M * M
        return 1.0 / self.n * 1.0 / self.n * \
                ( self.J / self.kT ) * \
                ( 1.0 / steps * np.sum( M2 ) - ( 1.0 / steps * np.sum( M ) ) * ( 1.0 / steps * np.sum( M ) ) )

if __name__ == '__main__':
    np.random.seed( 1 ) # use the same set of random numbers... forever!
    J, n = 1.5, 100
    nsteps = 20000
    heat = []
    magn = []
    susz = []
    tempRange = np.arange( 1.0, 4.0, 0.2 )
    for kT in tempRange:
        print "*"*50, kT
        lat = Lattice_old( J, n, kT )
        lat.get_initial_E_total() # I don't think this line actually does anything... 
        # To check, uncomment the following and try it.
#        print lat.get_initial_E_total(), "hey!!"        
        for step in range( nsteps ):
            lat.relax( step )
        heat.append( lat.heat( nsteps ) )
        #magn.append(np.sum(lat.magn(nsteps)))
        magn.append( np.sum( lat.magnet( nsteps ) ) )
        susz.append( lat.susz( nsteps ) )
    print heat
    
    # TODO: implement this with pyplot
    import pylab
    print "Simulation done."
    pylab.plot( tempRange, heat, 'ro' )
    pylab.savefig( 'heat.png' )
    print "Save heat done."
    pylab.clf()
    pylab.plot( tempRange, magn, 'ro' )
    pylab.savefig( 'magn.png' )
    print "Save magn done."
    pylab.clf()
    pylab.plot( tempRange, susz, 'ro' )
    pylab.savefig( 'susz.png' )
    print "Save susz done."
