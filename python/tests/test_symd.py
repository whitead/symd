import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import probplot
from scipy.stats import normaltest
import matplotlib as mpl
from math import ceil
from symd import Symd

mpl.use('Agg')


class Test:
    summary = ''

    def print_summary(self, force_print=False):
        if not force_print and self.passed:
            print(f"{self.name}:Passed ({self.summary})")
        else:
            print("%s:Failed" % self.name)
            print(self.summary)
            try:
                print('Check plots:')
                for p in self.plots:
                    print(f'\t{p.path}')
            except AttributeError:
                pass


class Plot:
    def __init__(self, name, path, caption="", size=[8, 6]):
        self.name = name
        self.path = path
        self.caption = caption if caption else name
        self.size = size


class Tester:
    def __init__(self, md, name):
        self.name = name
        self.tests = []
        self.md = md
        md.remove_overlap()
        md.shrink()
        self.add(ExecuteTest(md, name))
        if md.runParams['temperature'] == 0:
            md.log_output(period=1)
            self.add(EnergyConservationTest(md, name))
        else:
            md.log_output(period=100)
            self.add(ThermostatTest(md, name))

    def add(self, test):
        self.tests.append(test)

    def run(self):
        print("SimpleMD Tester `%s`:" % self.name)
        success = True
        for t in self.tests:
            success = success and t.run()
            print("*", end=" ")
            t.print_summary()
        if success:
            self.md.clean_files()


class EnergyConservationTest(Test):
    def __init__(self, md, name):
        self.short_name = "Energy_%s" % name
        self.name = "Energy Conservation `%s`" % name
        self.md = md
        self.passed = False

        self.dtol = 1e-10
        self.ptol = 0.05
        self.equil_frac = 0.3

    def run(self):

        if not md.executed:
            md.run()

        start_index = round(len(md.te) * self.equil_frac)
        energy = md.te[start_index:]
        if not np.all(np.isfinite(energy)):
            raise ValueError('NaN in energy')
        mean = np.mean(energy)
        dev = np.std(energy)

        # Check if it's changing a very small amount

        if dev < self.dtol:
            self.passed = True
            self.summary = "Energy fluctuation is %g" % dev
            return self.passed

        # Check if there's a correlation over time with energy (not good)

        norm_energy = (energy - mean) / dev
        cor = pearsonr(list(range(len(energy))), norm_energy)

        if cor[1] > self.ptol:
            self.passed = True

        # Make plots showing the energy change over time and compare its
        # distribution with normal

        self.summary = "Time Correlation = %g (p-value: %g)" % cor

        self.plots = [
            Plot(
                "%s Energy distribution plot" % self.name,
                "%s_energy_dist.png" % self.short_name,
            )
        ]

        fig = plt.figure(figsize=(12, 6))

        ax = plt.subplot2grid((2, 2), (0, 0), colspan=2)

        ax.plot(md.out_times[start_index:], energy)

        # take running average
        window = ceil(len(energy) / 10.0)
        weights = np.ones(window)
        energy_smooth = np.convolve(
            weights / weights.sum(), md.te[(start_index - window + 1):], mode="valid"
        )

        ax.plot(md.out_times[start_index:],
                energy_smooth, linewidth=2, color="r")
        ax.set_title("Time vs Energy")
        ax.set_xlabel("time")
        ax.set_ylabel("energy")
        ax.grid(True)

        ax = plt.subplot2grid((2, 2), (1, 0))
        ax.hist(norm_energy, 50, density=True,
                facecolor="yellow", alpha=0.75)

        ax.set_xlabel("Normalized Energy")
        ax.set_ylabel("Density")

        ax.grid(True)

        ax = plt.subplot2grid((2, 2), (1, 1))

        ((osm, osr), (slope, intercept, r)) = probplot(norm_energy, fit=True)

        ax.plot(osm, osr, "o", osm, slope * osm + intercept)
        ax.set_title("Probability Plot")
        ax.set_xlabel("Quantiles")
        ax.set_ylabel("Ordered Values")

        fig.tight_layout()

        plt.savefig(self.plots[0].path)

        return self.passed


class ExecuteTest(Test):
    def __init__(self, md, name):
        self.name = "Execute test `%s`" % name
        self.md = md
        self.passed = False

    def run(self):
        md.run()
        if md.outerr != b"":
            self.summary = "Execute test stderr out:\n %s" % md.outerr
        else:
            self.passed = True
        return self.passed


class ThermostatTest(Test):
    """
    Checks to see if the temperature is the same as the set temperature
    """

    def __init__(self, md, name):
        self.short_name = "thermtest_%s" % name
        self.name = "Thermostat test `%s`" % name
        self.md = md
        self.passed = False
        self.equil_frac = 0.30

        self.mtol = 1e-10
        self.ptol = 0.05

    def run(self):
        # check to see if there is a temperature
        if md.runParams["temperature"] == 0:
            self.passed = True
            self.summary = "No set temperature"
            return self.passed

        if not md.executed:
            md.run()

        start_index = round(len(md.temperature) * self.equil_frac)

        temperature = md.temperature[start_index:]
        set_temperature = md.runParams["temperature"]

        # Check if it's solid as a rock less than the set.
        if (np.mean(temperature) - set_temperature) / set_temperature < self.mtol:
            self.passed = True
            self.summary = "Temperature ppm difference less than %g" % self.mtol
            return self.passed

        # If it's not, let's check to see if the deviation is systematic and
        # see if the the deviation is normally distributed

        self.plots = [
            Plot(
                "%s Temperature distribution plot" % self.name,
                "%s_temp_dist.png" % self.short_name,
            )
        ]

        temperature_deviation = temperature - set_temperature
        k2, pvalue = normaltest(temperature_deviation)
        self.summary = "Deviation is %g %%. Deviation normality p-value = %g" % (
            (np.mean(temperature) - set_temperature) / set_temperature,
            pvalue,
        )

        # I've flipped this - it should be Maxwell-Boltzmann distributed not normal!
        if pvalue < self.ptol:
            self.passed = True

        fig = plt.figure(figsize=(12, 6))

        ax = plt.subplot2grid((2, 2), (0, 0), colspan=2)

        ax.plot(md.out_times[start_index:], temperature)

        # take running average
        window = ceil(len(temperature) / 10.0)
        weights = np.ones(window)
        temperature_smooth = np.convolve(
            weights / weights.sum(),
            md.temperature[max(start_index - window + 1, 0):],
            mode="valid",
        )

        ax.plot(
            md.out_times[(len(md.out_times) - len(temperature_smooth)):],
            temperature_smooth,
            linewidth=2,
            color="r",
        )

        ax.axhline(set_temperature, color="green", linewidth=2)

        ax.set_title("Time vs Temperature")
        ax.set_xlabel("time")
        ax.set_ylabel("$T$")
        ax.grid(True)

        ax = plt.subplot2grid((2, 2), (1, 0))

        ax.hist(
            temperature - set_temperature, 50, density=True, facecolor="blue", alpha=0.75
        )

        ax.axvline(0, color="green", linewidth=2)

        ax.set_xlabel("$\Delta$Temperature")
        ax.set_ylabel("Density")
        ax.grid(True)

        ax = plt.subplot2grid((2, 2), (1, 1))

        ((osm, osr), (slope, intercept, r)) = probplot(
            temperature - set_temperature, fit=True
        )

        ax.plot(osm, osr, "o", osm, slope * osm + intercept)
        ax.set_title("Probability Plot")
        ax.set_xlabel("Quantiles")
        ax.set_ylabel("Ordered Values")

        fig.tight_layout()
        plt.savefig(self.plots[0].path)


if __name__ == "__main__":

    np.random.seed(0)

    cell = [15, 15]
    for i in range(1, 3):
        md = Symd(nparticles=5, cell=cell, ndims=2, force='soft',
                  group=i, steps=500, exeDir="test")
        runner = Tester(md, f"SF_2_{i}_NVE")
        runner.run()

        md = Symd(nparticles=5, cell=cell, ndims=2,
                  group=i, steps=500, exeDir="test")
        runner = Tester(md, f"LJ_2_{i}_NVE")
        runner.run()

        md = Symd(nparticles=5, cell=cell, temperature=0.5, ndims=2,
                  group=i, steps=5000, exeDir="test")
        runner = Tester(md, f"LJ_2_{i}_NVT")
        runner.run()
