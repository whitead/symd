import os
import sys
import subprocess
import json
from math import *
import numpy as np
import build_lattice

SimpleMDLocation = "../build/"


class SimpleMD:
    """SimpleMD Engine set-up and output parser"""

    def __init__(
        self,
        nparticles,
        ndims,
        steps=100,
        temperature=None,
        exeDir="",
        exePrefix=SimpleMDLocation,
        thermostat="bussi",
        force="lj",
    ):

        self.ndims = ndims
        self.nparticles = nparticles

        self.prefix = exeDir
        if self.prefix != "" and (not os.path.exists(self.prefix)):
            os.makedirs(self.prefix)
            self.createdDir = True
        else:
            self.createdDir = False

        self.executed = False
        self.do_log_output = False

        self.runParams = {
            "steps": steps,
            "n_dims": self.ndims,
            "n_particles": nparticles,
            "time_step": 0.005,
            "velocity_seed": 435423,
            "anderson_nu": 100,
            "bussi_taut": 5,
            "thermostat_seed": 54344,
            "rcut": 3,
            "skin": 3 * 0.2,
            "com_remove_period": 1000,
            # this is the initial velocity temperature if NVE is chosen
            "temperature": 0,
            "thermostat": thermostat,
            "force_type": force,
            "final_positions": os.path.join(self.prefix, "final_positions.xyz"),
            "print_period": max(1, (steps / 100))
        }

        self.exePrefix = exePrefix
        self.exe = os.path.join(exePrefix, 'simple-md')
        if temperature != None:
            self.runParams["temperature"] = temperature

        if force == "lj":
            self.runParams["lj_epsilon"] = 1
            self.runParams["lj_sigma"] = 1

    def log_positions(self, filename="positions.xyz", period=0):
        """enable logging of the xyz positions of the simulation. Default is to output 100 frames"""
        if period == 0:
            period = ceil(self.runParams["steps"] / 100.0)

        self.runParams["position_log_period"] = int(period)
        self.runParams["positions_log_file"] = os.path.join(
            self.prefix, filename)

    def log_output(self, filename="md.log", period=0):
        if period != 0:
            self.runParams["print_period"] = period
        self.runParams["log_file"] = os.path.join(self.prefix, filename)
        self.do_log_output = True

    def setup_masses(self, masses=1, masses_file=None):
        """Creates a masses file. The masses variable is expanded to be the
        same size as the number of particles"""

        if masses_file == None:
            masses_file = 'masses.txt'
        self.runParams["masses_file"] = masses_file

        try:
            i = 0
            while len(masses) < self.nparticles:
                masses.append(masses[i])
                i += 1
        except TypeError:
            # an int was passed, not array
            masses = [masses for x in range(self.nparticles)]

        with open(self.runParams["masses_file"], "w") as f:
            for m in masses:
                f.write("%d\n" % m)

    def setup_positions(
        self,
        density=None,
        box_size=None,
        start_positions=None,
        overwrite=True,
        removeOverlap=False,
    ):
        """build a uniform lattice of particles given the box size or density"""

        if start_positions == None:
            start_positions = os.path.join(self.prefix, "start_positions.xyz")

        if density != None:
            volume = self.nparticles / density
            edge_size = volume ** (1.0 / self.ndims)
            box_size = [edge_size for x in range(self.ndims)]

        if len(box_size) != self.ndims:
            raise Exception(
                "Incorrect number of box dimensions. Must be %d" % self.ndims
            )

        self.runParams["start_positions"] = start_positions
        self.runParams["box_size"] = box_size

        if overwrite:
            increment = build_lattice.increment_size(
                self.ndims, box_size, self.nparticles
            )

            temp = sys.stdout

            with open(start_positions, "w") as f:
                sys.stdout = f

                build_lattice.enumerate_grid(
                    lambda x, y: build_lattice.print_grid(
                        x, y, increment, self.nparticles
                    ),
                    self.ndims - 1,
                    [int(ceil(x / increment)) for x in box_size],
                    [],
                )

            sys.stdout = temp

        if removeOverlap:
            # we create another MD engine to run with a soft potential and take the results

            overlapMD = SimpleMD(
                self.nparticles,
                self.ndims,
                temperature=None,
                exeDir=os.path.join(self.prefix, "overlap"),
                force="soft",
                exePrefix=self.exePrefix,
            )
            overlapMD.setup_positions(
                box_size=box_size,
                start_positions=start_positions,
                overwrite=False,
                removeOverlap=False,
            )
            overlapMD.runParams["steps"] = 1000
            overlapMD.runParams["print_period"] = 1000
            overlapMD.runParams["time_step"] = 0.001
            overlapMD.runParams["final_positions"] = os.path.join(
                self.prefix, "new_start_positions.xyz")
            if(not overlapMD.execute()):
                raise ValueError('Failed to create non-overlapping particles')
                # prevent overlap from cleaning out position file
            overlapMD.runParams["final_positions"] = ""
            overlapMD.clean_files()
            os.rename(os.path.join(
                self.prefix, "new_start_positions.xyz"), self.runParams['start_positions'])

        self.position_ready = True

    def clean_files(self):
        files = [
            "masses_file",
            "start_positions",
            "positions_log_file",
            "velocities_log_file",
            "forces_log_file",
            "log_file",
            "final_positions",
        ]
        for k in files:
            if k in self.runParams:
                if os.path.exists(self.runParams[k]):
                    os.remove(self.runParams[k])

        if self.createdDir:
            os.rmdir(self.prefix)

    def execute(self):

        if not self.position_ready:
            raise Exception("Must call setup_positions() before executing")

        output_lines = ceil(
            self.runParams["steps"] / self.runParams["print_period"])

        self.out_times = np.empty(output_lines)
        self.temperature = np.empty(output_lines)
        self.ke = np.empty_like(self.temperature)
        self.pe = np.empty_like(self.temperature)
        self.te = np.empty_like(self.temperature)
        self.htherm = np.empty_like(self.temperature)
        self.v = np.empty_like(self.temperature)

        out_arrays = (self.temperature, self.pe, self.ke, self.te, self.htherm, self.v)

        proc = subprocess.Popen(
            self.exe,
            stdin=subprocess.PIPE,
            bufsize=4096,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        output, outerr = proc.communicate(json.dumps(self.runParams).encode())
        if proc.returncode != 0:
            print(outerr.decode())
            raise RuntimeError(
                f'Failed to complete simulation. Ended with code {proc.returncode}')

        if self.do_log_output:
            with open(self.runParams["log_file"], "wb") as f:
                f.write(output)
        self.outerr = outerr
        index = 0

        for line in output.decode().split("\n"):
            try:
                sline = line.split()
                if len(sline) < 2:
                    continue
                step = int(sline[0])
                self.out_times[index] = float(sline[1])
                for i, oa in zip(list(range(len(out_arrays))), out_arrays):
                    # +1 for step and +1 for time
                    oa[index] = float(sline[i + 2])
                index += 1
            except ValueError:
                pass
        self.executed = True
        return True
