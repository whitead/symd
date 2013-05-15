import os
import build_lattice

SimpleMDLocation = "../"

class SimpleMD:
    """SimpleMD Engine set-up and output parser"""    

    def __init__(self, nparticles, ndims,
                 exeDir='',
                 exePrefix=SimpleMDLocation,
                 integrator='vverlet', 
                 thermostat='bussi', 
                 force='lj'):

        self.name = '%s_%s_%s' % (force, integrator, thermostat)

        self.ndims = ndims
        self.nparticles = nparticles

        self.prefix = exeDir
        if (self.prefix != "" and (not os.path.exists(self.prefix))):
            os.makedirs(self.prefix)

        self.runParams = { 'steps':100,
                           'nparticles':nparticles,
                           'time_step':0.005,
                           'velocity_seed': 435423,
                           'anderson_nu':100,
                           'bussi_taut':5,
                           'thermostat_seed':54344,
                           'rcut':3,
                           'print_period':1,
                           'masses_file':self.prefix + os.sep + 'masses.txt' }


    def setup_masses(masses = 1, masses_file=self.runParams['masses_file']):        
        """Creates a masses file. The masses variable is expanded to be the 
        same size as the number of particles"""

        self.runParams['masses_file'] = masses_file

        try:
            i = 0
            while(len(masses) < self.nparticles):
                masses.append(masses[i])
                i += 1
        except TypeError:
            #an int was passed, not array
            masses = [masses for x in range(self.nparticles)]

        with open(self.runParams['masses_file'], 'w') as f:
            for m in masses:
                f.writeline("%d\n" % m)
            
    def setup_positions(box_size, start_positions=None, overwrite=(start_positions == None)):
        """build a uniform lattice of particles given the box size or density"""

        if(start_positions == None):
            start_positions = self.prefix + os.sep + 'start_positions.xyz'

        if(len(box_size) != self.ndims):
            raise Exception('Incorrect number of box dimensions. Must be %d' % self.ndims)

        self.runParams['start_positions'] = start_positions

        for i in range(self.ndims):
            self.runParams['box_%d_size'] = box_size[i]

        if(overwrite):
            increment = increment_size(n_dims, box_size, number)

            temp = sys.stdout            

            with open(start_positions, 'w') as f:
                sys.stdout = f

                build_lattice.enumerate_grid(
                    lambda x,y: build_lattice.print_grid(x, y, increment, number), 
                    n_dims - 1,
                    [int(ceil(x / increment)) for x in box_size], 
                    [])

            sys.stdout = temp
            


md = new SimpleMD(100, 3, "test")
md.setup_masses()
