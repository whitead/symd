from SimpleMD import SimpleMD
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import probplot
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from math import ceil

class Test:

    def print_summary(self, force_print=False):
        if(not force_print and self.passed):
            print "%s:Passed" % self.name
        else:
            print "%s:Failed" % self.name
            print self.summary
            for p in self.plots:
                print "![%s](%s \"%s\")" % (p.caption, p.path, p.name)
                print p.caption
                print "\n"


class Plot:
    def __init__(self, name, path, caption="", size=[8,6]):
        self.name = name
        self.path = path
        self.caption = caption
        self.size = size


class SanityTester:

    def __init__(self, md, name):
        self.tests = []
        self.add(EnergyConservationTest(md, name))
        self.add(EnergyConservationTest(md, name))

    def add(self, test):
        self.tests.append(test)

    def run(self):
        print "Sanity Tester:"
        for t in self.tests:
            success = t.run()
            print "*",
            t.print_summary()
            if(not success):
                return

class EnergyConservationTest(Test):

    def __init__(self, md, name):
        self.short_name = "Energy_%s" % name
        self.name = "Energy Conservation [%s]" % name
        self.md = md
        self.passed = False

        self.dtol = 1e-10
        self.ptol = 0.05
        self.equil_frac = 0.3

    def run(self):

        if(not md.executed):
            md.execute()
            #md.clean_files()
        
        start_index = round(len(md.te) * self.equil_frac)
        energy = md.te[start_index:]
        mean = np.mean(energy)
        dev = np.std(energy)

        if(dev < self.dtol):
            self.passed = True
            self.summary = "Energy fluctation is %g" % dev
            return self.passed

        norm_energy = (energy - mean) / dev 
        cor = pearsonr(range(len(energy)), norm_energy)

        if(abs(cor[1]) > self.ptol):
            self.passed = True

        self.summary = "Time Correlation = %g (p-value: %g)" % cor

        self.plots = [Plot("%s Q-Q Norm Plot" % self.name, "%s_energy_qq.png" % self.short_name)]

        fig = plt.figure(figsize=(12, 6))

        ax = plt.subplot2grid((2,2), (0,0), colspan=2)

        ax.plot(md.out_times[start_index:], energy)
        
        #take running average
        window = ceil(len(energy) / 10.)
        weights = np.ones(window, 'd')
        energy_smooth = np.convolve(weights / weights.sum(), md.te[(start_index - window + 1):], mode='valid')

        ax.plot(md.out_times[start_index:], energy_smooth, linewidth=2, color="r")
        ax.set_title('Time vs Energy')
        ax.set_xlabel('time')
        ax.set_ylabel('energy')
        ax.grid(True)

        ax = plt.subplot2grid((2,2), (1,0))
        n, bins, patches = ax.hist(norm_energy,  50, normed=1, facecolor='yellow', alpha=0.75)

        ax.set_xlabel("Normalized Energy")
        ax.set_ylabel("Density")

        ax.grid(True)

        ax = plt.subplot2grid((2,2), (1,1))

        ((osm, osr), (slope, intercept, r)) = probplot(norm_energy, fit=True)

        ax.plot(osm, osr, 'o', osm, slope*osm + intercept)
        ax.set_title('Probability Plot')
        ax.set_xlabel('Quantiles')
        ax.set_ylabel('Ordered Values')

        xmin = np.amin(osm)
        xmax = np.amax(osm)
        ymin = np.amin(norm_energy)
        ymax = np.amax(norm_energy)
        posx = xmin + 0.70 * (xmax - xmin)
        posy = ymin + 0.01 * (ymax - ymin)

        fig.tight_layout()

        plt.savefig(self.plots[0].path)


        return self.passed

                  

md = SimpleMD(100, 3, steps=50000, temperature=0.7, exeDir="test")
md.log_output(period=10)
md.setup_positions(0.7)
md.setup_masses(1)

runner = SanityTester(md, "LJ_NVT")
runner.run()

md = SimpleMD(2, 3, steps=50000, exeDir="test")
md.log_output(period=10)
md.setup_positions(box_size=[4,4,4])
md.setup_masses(1)

runner = SanityTester(md, "LJ_NVE")
runner.run()
