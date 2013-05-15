import os, sys, subprocess
from math import *
import numpy as np
import build_lattice

SimpleMDLocation = ".."


class SimpleMD:
    """SimpleMD Engine set-up and output parser"""    

    def __init__(self, nparticles, ndims, temperature=None,
                 exeDir='',
                 exePrefix=SimpleMDLocation,
                 integrator='vverlet', 
                 thermostat='bussi', 
                 force='lj'):


        self.ndims = ndims
        self.nparticles = nparticles

        self.prefix = exeDir
        if (self.prefix != "" and (not os.path.exists(self.prefix))):
            os.makedirs(self.prefix)
            self.createdDir = True
        else:
            self.createdDir = False

        self.executed = False

        self.runParams = { 'steps':100,
                           'n_dims': self.ndims,
                           'n_particles':nparticles,
                           'time_step':0.005,
                           'velocity_seed': 435423,
                           'anderson_nu':100,
                           'bussi_taut':5,
                           'thermostat_seed':54344,
                           'rcut':3,
                           'skin':3 * 0.2,
                           'com_remove_period':1000,
                           #this is the initial velocity temperature if NVE is chosen
                           'temperature':1, 
                           'print_period':1,
                           'masses_file':self.prefix + os.sep + 'masses.txt' }
        
        if(temperature != None):
            self.runParams['temperature'] = temperature
            self.exe = '%s%s%s_%s_%s' % (exePrefix, os.sep, force, integrator, thermostat)
        else:
            self.exe = '%s%s%s_%s' % (exePrefix, os.sep, force, integrator)

        if(force == 'lj'):
            self.runParams['lj_epsilon'] = 1
            self.runParams['lj_sigma'] = 1

    def setup_masses(self, masses = 1, masses_file = None):        
        """Creates a masses file. The masses variable is expanded to be the 
        same size as the number of particles"""
        
        if(masses_file == None):
            masses_file = self.runParams['masses_file']
        else:
            self.runParams['masses_file'] = self.prefix + os.sep + masses_file

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
                f.write("%d\n" % m)
            
    def setup_positions(self, box_size, start_positions=None, overwrite=True):
        """build a uniform lattice of particles given the box size or density"""

        if(start_positions == None):
            start_positions = self.prefix + os.sep + 'start_positions.xyz'
            

        if(len(box_size) != self.ndims):
            raise Exception('Incorrect number of box dimensions. Must be %d' % self.ndims)

        self.runParams['start_positions'] = start_positions

        for i in range(self.ndims):
            self.runParams['box_%d_size' % (i + 1)] = box_size[i]

        if(overwrite):
            increment = build_lattice.increment_size(self.ndims, box_size, self.nparticles)

            temp = sys.stdout            

            with open(start_positions, 'w') as f:
                sys.stdout = f

                build_lattice.enumerate_grid(
                    lambda x,y: build_lattice.print_grid(x, y, increment, self.nparticles), 
                    self.ndims - 1,
                    [int(ceil(x / increment)) for x in box_size], 
                    [])

            sys.stdout = temp

    def clean_files(self):
        files = ['masses_file', 'start_positions', 'positions_log_file',
                 'velocities_log_file', 'forces_log_file']
        for k in files:
            if(self.runParams.has_key(k)):
                if(os.path.exists(self.runParams[k])):
                    os.remove(self.runParams[k])

        if(self.createdDir):
            os.rmdir(self.prefix)
            
    def execute(self):
        output_lines = self.runParams['steps']

        self.temperature = np.empty(output_lines)
        self.ke = np.empty_like(self.temperature)
        self.pe = np.empty_like(self.temperature)
        self.te = np.empty_like(self.temperature)
        self.htherm = np.empty_like(self.temperature)

        out_arrays = (self.temperature, self.pe, self.ke, self.te, self.htherm)

        input_string = []
        
        for k in self.runParams.iterkeys():
            input_string.append('%s %s\n' % (k, self.runParams[k]))

        proc = subprocess.Popen(self.exe, stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE)
        
        output = proc.communicate(''.join(input_string))[0]

        for line in output.split('\n'):
            try:
                sline = line.split()
                step = int(sline[0])
                for i, oa in zip(range(len(out_arrays)), out_arrays):
                    oa[step] = float(sline[i + 2])
            except:
                pass

        self.executed = True
        


