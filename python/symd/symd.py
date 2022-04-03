import os
import sys
import subprocess
import json
from math import *
from weakref import WeakKeyDictionary
import numpy as np
from .groups import *

SymdLocation = ""


class Symd:
    """Symd Engine set-up and output parser"""

    def __init__(
        self,
        nparticles,
        ndims,
        cell=None,
        density=None,
        group=1,
        wyckoffs=None,
        images=None,
        steps=100,
        temperature=None,
        start_temperature=0.5,
        pressure=None,
        exeDir="",
        exePrefix=SymdLocation,
        thermostat="baoab",
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
        self.positions = None

        self.runParams = {
            "steps": steps,
            "n_particles": nparticles,
            "time_step": 0.005,
            "seed": 435423,
            "anderson_nu": 100,
            "langevin_gamma": 0.1,
            "bussi_taut": 0.5,
            "thermostat_seed": 54344,
            "rcut": 10,
            "skin": 3 * 0.2,
            "temperature": 0,
            "pressure": 0,
            "box_update_period": 0,
            "start_temperature": start_temperature,
            "thermostat": thermostat,
            "force_type": force,
            "final_positions": os.path.join(self.prefix, "final_positions.dat"),
            "cell_log_file": os.path.join(self.prefix, "cell_log_file.dat"),
            "print_period": max(1, (steps / 100)),
            "cell": cell,
        }

        self.exePrefix = exePrefix
        self.exe = os.path.join(exePrefix, f"symd{self.ndims}")
        if temperature != None:
            self.runParams["temperature"] = temperature
            self.runParams["start_temperature"] = temperature
        else:
            self.runParams["thermostat"] = None
        if pressure != None:
            self.runParams["pressure"] = pressure
            self.runParams["box_update_period"] = 5
        else:
            self.runParams["pressure"] = None
        if force == "lj":
            self.runParams["lj_epsilon"] = 1
            self.runParams["lj_sigma"] = 1

        if images:
            if type(images) == int:
                self.runParams["n_images"] = images
            elif type(images) == list:
                self.runParams["images"] = images

        if type(group) == int:
            gname = "group-" + str(group)
            group = load_group(group, ndims)
        else:
            gname = "group-custom"
        outputs = prepare_input(group, ndims, nparticles, gname, self.prefix)

        if wyckoffs is not None:
            if type(wyckoffs) == int:
                wyckoffs = [wyckoffs] * len(outputs)
            ws = []
            for i, w in enumerate(wyckoffs):
                if i >= len(outputs):
                    raise ValueError("Too many Wyckoff positions specified")
                ws.append({"group": outputs[i], "n_particles": w})
            self.runParams["wyckoffs"] = ws
        else:
            wyckoffs = []

        if cell is None:
            if density is None:
                raise ValueError("Must specify density or cell")
            else:
                cell = get_cell(density, group, ndims, nparticles, wyckoffs)
        else:
            if density is not None:
                raise ValueError("Cannot specify both density and cell")
        self.runParams["cell"] = cell
        self.runParams["start_positions"] = os.path.join(self.prefix, gname + ".dat")
        self.runParams["group"] = os.path.join(self.prefix, f"{gname}.json")
        # compute cell particle count
        self.cell_nparticles = cell_nparticles(group, nparticles, *wyckoffs)

    def remove_overlap(self):
        """Attempt to remove overlapping particles"""
        cache_params = self.runParams.copy()
        self.runParams["force_type"] = "soft"
        self.runParams["thermostat"] = "baoab"
        self.runParams["pressure"] = None
        self.runParams["box_update_period"] = 0
        self.runParams["temperature"] = 1.0
        self.runParams["steps"] = 1000
        self.runParams["final_positions"] = self.runParams["start_positions"]
        self.run()
        self.runParams["force_type"] = "lj"
        self.runParams["thermostat"] = "baoab"
        self.runParams["langevin_gamma"] = 10
        self.runParams["temperature"] = 0.01
        self.runParams["start_temperature"] = 0.5
        self.runParams["lj_sigma"] = 1.0
        self.runParams["steps"] = 1000
        self.runParams["final_positions"] = self.runParams["start_positions"]
        self.run()
        self.runParams = cache_params

    def shrink(self):
        """Run a brief NPT simulation"""
        cache_params = self.runParams.copy()
        self.runParams["force_type"] = "lj"
        self.runParams["thermostat"] = "baoab"
        self.runParams["langevin_gamma"] = 10
        self.runParams["temperature"] = 0.1
        self.runParams["start_temperature"] = 0.2
        self.runParams["steps"] = 5000
        self.runParams["final_positions"] = self.runParams["start_positions"]
        self.runParams["pressure"] = 0.1
        self.runParams["box_update_period"] = 3
        self.run()
        self.runParams = cache_params
        self.runParams["cell"] = self.read_cell()

    def log_positions(self, filename="positions.xyz", period=0, frames=0):
        """enable logging of the xyz positions of the simulation. Default is to output 100 frames"""
        if period == 0:
            if frames == 0:
                frames = 100
            period = ceil(self.runParams["steps"] / frames)

        self.runParams["position_log_period"] = int(period)
        self.runParams["positions_log_file"] = os.path.join(self.prefix, filename)

    def log_output(self, filename="md.log", period=0, frames=0):
        """Log simulation output"""
        if period == 0:
            if frames == 0:
                frames = 100
            period = ceil(self.runParams["steps"] / frames)
        self.runParams["print_period"] = int(period)
        self.runParams["log_file"] = os.path.join(self.prefix, filename)
        self.do_log_output = True

    def setup_masses(self, masses=1, masses_file=None):
        """Creates a masses file. The masses variable is expanded to be the
        same size as the number of particles"""

        if masses_file == None:
            masses_file = "masses.txt"
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

    def clean_files(self):
        """Remove output files"""
        files = [
            "masses_file",
            "start_positions",
            "positions_log_file",
            "velocities_log_file",
            "forces_log_file",
            "log_file",
            "final_positions",
            "cell_log_file",
            "group",
        ]
        for k in files:
            if k in self.runParams:
                if os.path.exists(self.runParams[k]):
                    os.remove(self.runParams[k])
        if "wyckoffs" in self.runParams:
            for w in self.runParams["wyckoffs"]:
                if os.path.exists(w["group"]):
                    os.remove(w["group"])
        if self.createdDir:
            os.rmdir(self.prefix)

    def run(self):
        """Run the simulation"""
        output_lines = ceil(self.runParams["steps"] / self.runParams["print_period"])

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
        if "positions_log_file" in self.runParams:
            try:
                self.read_positions()
            except AssertionError as e:
                print(e)
                print("could not read positions")

        # update so next run is valid
        self.runParams["cell"] = self.read_cell()
        self.runParams["start_positions"] = self.runParams["final_positions"]

        # Do it at the end in case we want to catch
        # and still use results
        if proc.returncode != 0:
            print(json.dumps(self.runParams))
            print(output.decode())
            print(outerr.decode())
            raise RuntimeError(
                f"Failed to complete simulation. Ended with code {proc.returncode}"
            )
        self.executed = True

        return True

    def read_cell(self, bravais=False):
        """Return the output unit cell.
        :param bravais: if True, return the cell fit to be in Bravais (projected) for use in p1 or non-symmetric MD. If False, return a unit cell that will repeat the simulation.
        """
        cell = []
        with open(self.runParams["cell_log_file"], "r") as f:
            for line in f.readlines():
                cell.extend([float(s) for s in line.split()])
        N = len(cell)
        if bravais:
            return cell[: N // 2]
        return cell[N // 2 :]

    def number_density(self):
        """returns the packing fraction of the system"""
        cell = self.read_cell(True)
        cell = np.array(cell).reshape((self.ndims, self.ndims))
        return self.cell_nparticles / cell_volume(cell)

    def read_positions(self):
        """Reads the positions from the positions log file"""
        # read one frame to get position conut
        with open(self.runParams["positions_log_file"], "r") as f:
            N = -1
            for line in f.readlines():
                if "Frame" in line:
                    if N > 0:
                        break
                    else:
                        N = 0
                else:
                    N += 1
        N -= 1
        with open(self.runParams["positions_log_file"], "r") as f:
            lines = f.readlines()
            K = len(lines)
            T = K // (N + 2)  # 2 for the header
            assert (N + 2) * T == K
            positions = np.empty((T, N, self.ndims))
            j = -1  # since we have Frame appear before first frame
            k = 0
            for i, line in enumerate(lines):
                sline = line.split()
                if len(sline) == 4:
                    positions[j, k, :] = [float(x) for x in sline[1 : self.ndims + 1]]
                    k += 1
                elif "Frame" in line:
                    j += 1
                    k = 0
        if self.positions is None:
            self.positions = positions
        else:
            self.positions = np.concatenate((self.positions, positions))
