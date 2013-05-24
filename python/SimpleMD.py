import os, sys, subprocess
from math import *
import numpy as np
import build_lattice
import gofr

SimpleMDLocation = ".."


class SimpleMD:
    """SimpleMD Engine set-up and output parser"""    

    def __init__(self, nparticles, ndims, steps=100, temperature=None,
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
        self.do_log_output = False

        self.runParams = { 'steps':steps,
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
                           'temperature':0, 
                           'print_period':max(1, (steps / 100)),
                           'masses_file':self.prefix + os.sep + 'masses.txt' }

        self.exePrefix = exePrefix
        if(temperature != None):
            self.runParams['temperature'] = temperature
            self.exe = '%s%s%s_%s_%s' % (exePrefix, os.sep, force, integrator, thermostat)
        else:
            self.exe = '%s%s%s_%s' % (exePrefix, os.sep, force, integrator)

        if(force == 'lj'):
            self.runParams['lj_epsilon'] = 1
            self.runParams['lj_sigma'] = 1

    def log_positions(self, filename='positions.xyz', period=0):
        """enable logging of the xyz positions of the simulation. Default is to output 100 frames"""
        if(period == 0):
            period = ceil(self.runParams['steps'] / 100.)

        self.runParams['position_log_period'] = int(period)
        self.runParams['positions_log_file'] = self.prefix + os.sep + filename


    def log_output(self, filename='md.log', period=0):
        if(period != 0):
            self.runParams['print_period'] = period
        self.runParams['log_file'] = self.prefix + os.sep + filename 
        self.do_log_output =True

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

        self.mass_ready = True
            
    def setup_positions(self, density=None, box_size=None, start_positions=None, overwrite=True, removeOverlap=True):
        """build a uniform lattice of particles given the box size or density"""

        if(start_positions == None):
            start_positions = self.prefix + os.sep + 'start_positions.xyz'
            
        if(density != None):
            volume = self.nparticles / density
            edge_size = volume ** (1. / self.ndims)
            box_size = [edge_size for x in range(self.ndims)]
        
        if(len(box_size) != self.ndims):
            raise Exception('Incorrect number of box dimensions. Must be %d' % self.ndims)

        self.runParams['start_positions'] = start_positions

        for i in range(self.ndims):
            self.runParams['box_%d_size' % (i + 1)] = box_size[i]

        self.box_size = box_size

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

        if(removeOverlap):
            #we create another MD engine to run with a soft potential and take the results

            overlapMD = SimpleMD(self.nparticles, self.ndims, temperature=None, 
                                 exeDir=self.prefix + os.sep + "overlap",
                                 force='soft',
                                 exePrefix=self.exePrefix)
            overlapMD.setup_positions(box_size=box_size, start_positions=start_positions, overwrite=False, removeOverlap=False)
            overlapMD.runParams['steps'] = 1000
            overlapMD.runParams['print_period'] = 1000
            overlapMD.runParams['time_step'] = 0.001
            overlapMD.runParams['positions_log_file'] = overlapMD.prefix + os.sep + 'positions.xyz'
            overlapMD.runParams['position_log_period'] = overlapMD.runParams['steps']
            overlapMD.setup_masses(1)
            overlapMD.execute()

            #get the output. Skip first 2 lines and atomic symbol
            with open(overlapMD.runParams['positions_log_file'], 'r') as infile:
                with open(start_positions, 'w') as outfile:
                    lines = infile.readlines()
                    for line in lines[2:]:
                        outfile.write("".join([x + " " for x in line.split()[1:(self.ndims + 1)]]))
                        outfile.write("\n")

            #prevent overlap from cleaning out position file
            overlapMD.runParams['start_positions'] = ''
            overlapMD.clean_files()

        self.position_ready = True
                    
            

    def clean_files(self):
        files = ['masses_file', 'start_positions', 'positions_log_file',
                 'velocities_log_file', 'forces_log_file', 'log_file']
        for k in files:
            if(self.runParams.has_key(k)):
                if(os.path.exists(self.runParams[k])):
                    os.remove(self.runParams[k])

        if(self.createdDir):
            os.rmdir(self.prefix)
            
    def execute(self):

        if(not self.mass_ready):
            raise Exception("Must call setup_masses() before executing")

        if(not self.position_ready):
            raise Exception("Must call setup_positions() before executing")


        output_lines = ceil(self.runParams['steps'] / self.runParams['print_period'])

        self.out_times = np.empty(output_lines)
        self.temperature = np.empty(output_lines)
        self.ke = np.empty_like(self.temperature)
        self.pe = np.empty_like(self.temperature)
        self.te = np.empty_like(self.temperature)
        self.htherm = np.empty_like(self.temperature)

        out_arrays = (self.temperature, self.pe, self.ke, self.te, self.htherm)

        input_string = []
        
        for k in self.runParams.iterkeys():
            input_string.append('%s %s\n' % (k, self.runParams[k]))

        proc = subprocess.Popen(self.exe, stdin=subprocess.PIPE, bufsize=4096,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        output,outerr = proc.communicate(''.join(input_string))

        if(self.do_log_output):
            with open(self.runParams['log_file'], 'w') as f:
                f.write(output)

        self.outerr = outerr

        index = 0

        for line in output.split('\n'):
            try:
                sline = line.split()
                if(len(sline) < 2):
                    continue
                step = int(sline[0])
                self.out_times[index] = float(sline[1])
                for i, oa in zip(range(len(out_arrays)), out_arrays):
                    # +1 for step and +1 for time
                    oa[index] = float(sline[i + 2]) 
                index += 1
            except ValueError:
                pass

        self.executed = True

    def calc_gofr(self, bin_width=0.1):
        """Calculates g(r)""" 
        #check that we have structures and a cube
        if(not self.executed):
            raise Exception("Must execute SimpleMD before calculating statistics")
        if(not self.runParams.has_key('positions_log_file')):
            raise Exception("Must call log_positions with SimpleMD to calculate gofr")
        return gofr.calc_gofr(self.runParams['positions_log_file'], self.box_size, 
                         bin_width=bin_width, ndims=self.ndims)
        


